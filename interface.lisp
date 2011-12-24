(in-package "NAPA-FFT")

#+sb-thread
(progn
  (declaim (type sb-thread:mutex *factors-lock*))
  (defvar *factors-lock* (sb-thread:make-mutex)))

(declaim (type complex-sample-array *fwd-factors* *bwd-factors*))
(defvar *fwd-factors* (make-all-factors 20 1))
(defvar *bwd-factors* (make-all-factors 20 -1))

(defun %ensure-factor-size (symbol size direction)
  (assert (= 1 (logcount size)))
  (let ((factors (symbol-value symbol)))
    (when (>= (length factors) (* 2 size))
      (return-from %ensure-factor-size factors)))
  (flet ((update ()
           (let ((factors (symbol-value symbol)))
             (cond ((>= (length factors) (* 2 size))
                    factors)
                   (t
                    (setf (symbol-value symbol)
                          (make-all-factors (integer-length (1- size))
                                            direction
                                            factors)))))))
    #+sb-thread
    (sb-thread:with-mutex (*factors-lock*)
      (update))
    #-sb-thread (update)))

(defmacro ensure-factor-size (symbol size direction)
  (let ((_factors   (gensym "FACTORS"))
        (_size      (gensym "SIZE"))
        (_direction (gensym "DIRECTION")))
    `(let ((,_factors ,symbol)
           (,_size    ,size)
           (,_direction ,direction))
       (declare (type complex-sample-array ,_factors)
                (type index ,_size)
                (type (member -1 1) ,_direction))
       (if (>= (truncate (length ,_factors) 2)
               size)
           ,_factors
           (%ensure-factor-size ',symbol ,_size ,_direction)))))

(defvar *forward*)
(defvar *forward-scale*)
(defvar *backward*)
(defvar *backward-scale*)

(macrolet ((fill-vectors (lb)
             (let ((max (integer-length array-total-size-limit)))
               `(progn
                  (declaim (type (simple-array fft-function (,(1+ max)))
                                 *forward* *backward*))
                  (declaim (type (simple-array (fft-function :scale t) (,(1+ max)))
                                 *forward-scale* *backward-scale*))
                  (multiple-value-bind (fwd fwd-scale bwd bwd-scale)
                    (build-fft-routine-vectors ,lb ,max)
                    (setf *forward*        fwd
                          *forward-scale*  fwd-scale
                          *backward*       bwd
                          *backward-scale* bwd-scale))
                  (values)))))
  ;; specializing FFTs up to 2^8 seem to be plenty
  (fill-vectors 8))

(defun find-fft-function (size direction &key (scale nil))
  (check-type size (and index (satisfies power-of-two-p)))
  (check-type direction (member :forward :backward))
  (check-type scale (or null double-float))
  (multiple-value-bind (functions factors)
      (ecase direction
        (:forward
         (values (if scale
                     *forward-scale*
                     *forward*)
                 (ensure-factor-size *fwd-factors* size 1)))
        (:backward
         (values (if scale
                     *backward-scale*
                     *backward*)
                 (ensure-factor-size *bwd-factors* size -1))))
    (let ((lb (integer-length (1- size))))
      (values (aref functions lb) factors))))

(defun %fft (size dst dst-offset src src-offset direction
             tmp1 tmp2 tmp-offset)
  (declare (type index size src-offset dst-offset tmp-offset)
           (type complex-sample-array dst src tmp1 tmp2))
  (multiple-value-bind (vector factors)
      (ecase direction
        (1
         (values *forward*  (ensure-factor-size *fwd-factors* size 1)))
        (-1
         (values *backward* (ensure-factor-size *bwd-factors* size -1))))
    (let ((len (integer-length (1- size))))
      (assert (= size (ash 1 len)))
      (funcall (aref vector len)
               factors size
               dst     dst-offset
               src     src-offset
               tmp1 tmp2 tmp-offset))))

(defun %fft-scale (size dst dst-offset src src-offset direction
                   scale
                   tmp1 tmp2 tmp-offset)
  (declare (type index size src-offset dst-offset tmp-offset)
           (type double-float scale)
           (type complex-sample-array dst src tmp1 tmp2))
  (multiple-value-bind (vector factors)
      (ecase direction
        (1
         (values *forward-scale*  (ensure-factor-size *fwd-factors* size 1)))
        (-1
         (values *backward-scale* (ensure-factor-size *bwd-factors* size -1))))
    (let ((len (integer-length (1- size))))
      (assert (= size (ash 1 len)))
      (funcall (aref vector len)
               factors size
               dst     dst-offset
               src     src-offset
               scale
               tmp1 tmp2 tmp-offset))))

(defstruct (fft-instance
            (:constructor %make-fft-instance))
  (size        0 :type index)
  (direction nil :type (member :forward :backward))
  (function  nil :type (or fft-function (fft-function :scale t)))
  (factors   nil :type complex-sample-array)
  (scale     1d0 :type double-float)
  (tmp1      nil :type complex-sample-array)
  (tmp2      nil :type complex-sample-array))

(defun reset-fft-instance (instance &key size direction (scale nil scale-p))
  (check-type size (or null
                       (and index
                            (satisfies power-of-two-p))))
  (check-type direction (or null
                            (member :forward :backward)))
  (check-type scale (or null double-float))
  (if size
      (setf (fft-instance-size instance) size)
      (setf size (fft-instance-size instance)))
  (if direction
      (setf (fft-instance-direction instance) direction)
      (setf direction (fft-instance-direction instance)))
  (if (not scale-p)
      (setf scale (fft-instance-scale instance))
      (setf (fft-instance-scale instance)
            (or scale
                (setf scale
                      (ecase direction
                        (:forward 1d0)
                        (:backward (/ (float size 1d0))))))))
  (multiple-value-bind (function factors)
      (find-fft-function size direction :scale (/= scale 1d0))
    (let ((tmp1 (fft-instance-tmp1 instance))
          (tmp2 (fft-instance-tmp2 instance)))
      (when (< (length tmp1) size)
        (setf tmp1 (make-array size :element-type 'complex-sample)))
      (when (< (length tmp2) size)
        (setf tmp2 (make-array size :element-type 'complex-sample)))
      (setf (fft-instance-size      instance) size
            (fft-instance-direction instance) direction
            (fft-instance-function  instance) function
            (fft-instance-factors   instance) factors
            (fft-instance-scale     instance) scale
            (fft-instance-tmp1      instance) tmp1
            (fft-instance-tmp2      instance) tmp2)
      instance)))

(defun make-fft-instance (size direction &optional scale)
  (check-type size (and index
                        (satisfies power-of-two-p)))
  (check-type direction (member :forward :backward))
  (check-type scale (or null double-float))
  (unless scale
    (setf scale
          (ecase direction
            (:forward 1d0)
            (:backward (/ (float size 1d0))))))
  (multiple-value-bind (function factors)
      (find-fft-function size direction :scale (/= scale 1d0))
    (%make-fft-instance :size     size
                        :direction direction
                        :function function
                        :factors  factors
                        :scale    scale
                        :tmp1     (make-array size :element-type 'complex-sample)
                        :tmp2     (make-array size :element-type 'complex-sample))))

(defun fft (instance src &optional dst)
  (check-type instance fft-instance)
  (check-type src      complex-sample-array)
  (check-type dst      (or null complex-sample-array))
  (assert (= (length src) (fft-instance-size instance)))
  (unless dst
    (setf dst (make-array (length src) :element-type 'complex-sample)))
  (assert (= (length dst) (fft-instance-size instance)))
  (if (eql 1d0 (fft-instance-scale instance))
      (funcall (fft-instance-function instance)
               (fft-instance-factors instance) (fft-instance-size instance)
               dst 0
               src 0
               (fft-instance-tmp1 instance) (fft-instance-tmp2 instance) 0)
      (funcall (fft-instance-function instance)
               (fft-instance-factors instance) (fft-instance-size instance)
               dst 0
               src 0
               (fft-instance-scale instance)
               (fft-instance-tmp1 instance) (fft-instance-tmp2 instance) 0)))

(declaim (type (or null fft-instance) *fft-instance*))
(defvar *fft-instance* nil)

(defun sfft (src &optional direction dst)
  (unless direction
    (setf direction :forward))
  (unless (typep src 'complex-sample-array)
    (setf src
          (typecase src
            ((simple-array double-float 1)
             (map 'complex-sample-array #'complex src))
            ((simple-array single-float 1)
             (map 'complex-sample-array (lambda (x)
                                          (complex (float x 1d0)))
                  src))
            (t
             (map 'complex-sample-array
                  (lambda (x)
                    (coerce x '(complex double-float)))
                  src)))))
  (check-type direction (member :forward :backward))
  (check-type dst (or null complex-sample-array))
  (let* ((size     (length src))
         (instance (setf *fft-instance*
                         (if *fft-instance*
                             (reset-fft-instance
                              *fft-instance*
                              :size size
                              :direction direction
                              :scale nil)
                             (make-fft-instance size direction)))))
    (fft instance src dst)))
