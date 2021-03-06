(in-package "NAPA-FFT")

(defmacro define-inline-function (name (&rest args) &body body)
  `(progn
     (declaim (inline ,name))
     (defun ,name (,@args)
       ,@body)
     (define-compiler-macro ,name (&rest arg-forms)
       `((lambda (,@',args) ,@',body)
         ,@arg-forms))))

(define-inline-function mul+i (x)
  (declare (type complex-sample x))
  #+ (and sbcl complex-float-vops)
  (sb-vm::swap-complex (conjugate x))
  #- (and sbcl complex-float-vops)
  (* x #c(0 1d0)))

(define-inline-function mul-i (x)
  (declare (type complex-sample x))
  #+ (and sbcl complex-float-vops)
  (conjugate (sb-vm::swap-complex x))
  #- (and sbcl complex-float-vops)
  (* x #c(0 -1d0)))

(define-inline-function mul+/-sqrt+i (x scale)
  (declare (type complex-sample x)
           (type double-float scale))
  #+ (and sbcl complex-float-vops)
  (let ((x (* x scale)))
    (+ x (sb-vm::swap-complex (conjugate x))))
  #- (and sbcl complex-float-vops)
  (* x (complex scale scale)))

(define-inline-function mul+/-sqrt-i (x scale)
  (declare (type complex-sample x)
           (type double-float scale))
  #+ (and sbcl complex-float-vops)
  (let ((x (* x scale)))
    (- x (sb-vm::swap-complex (conjugate x))))
  #- (and sbcl complex-float-vops)
  (* x (complex scale (- scale))))

(defun make-cooley-tukey-factors (size1 size2 direction
                                   &optional (coeffs (make-array
                                                      (* size1 size2)
                                                      :element-type 'complex-sample))
                                     (offset 0))
  (declare (type half-index size1 size2)
           (type (member -1 1) direction)
           (type complex-sample-array coeffs)
           (type index offset))
  (let* ((size   (* size1 size2))
         (base   (/ (* direction -2 pi) size)))
    (declare (type double-float base)
             (optimize speed))
    (dotimes (i size2)
      (let ((base (* base i)))
        (dotimes (j size1)
          (let ((k (+ (* j size2) i))
                (factor (* base j)))
            (setf (aref coeffs (+ offset k)) (cis factor))))))
    coeffs))

(declaim (type (integer -1) +factor-bias+))
(defconstant +factor-bias+ -1)

(defun make-all-factors (log-max-size direction &optional previous)
  (declare (type (or null complex-sample-array) previous))
  (when previous
    (let ((length (- (length previous) +factor-bias+)))
      (assert (power-of-two-p length))
      (when (>= (integer-length (1- length)) (1+ log-max-size))
        (return-from make-all-factors previous))))
  (let ((all-factors (make-array (+ (* 2 (ash 1 log-max-size))
                                    +factor-bias+)
                                 :element-type 'complex-sample)))
    (cond (previous
           (replace all-factors previous))
          (t
           (fill all-factors #c(0d0 0d0) :end (1+ +factor-bias+))))
    (loop for i from (if previous
                         (integer-length (1- (length previous)))
                         0)
            upto log-max-size do
      (let* ((total-size (ash 1 i))
             (size1      (ash 1 (truncate i 2)))
             (size2      (ash 1 (truncate (1+ i) 2))))
        (make-cooley-tukey-factors size1 size2 direction all-factors
                                   (+ total-size +factor-bias+))))
    all-factors))

(defun one-point-fft (twiddle size dst dst-offset src src-offset
                      tmp1 tmp2 tmp-offset)
  (declare (ignore twiddle size tmp1 tmp2 tmp-offset)
           (type complex-sample-array dst src)
           (type index dst-offset src-offset)
           (optimize speed (safety 0)))
  (setf (aref dst dst-offset) (aref src src-offset))
  dst)

(defun one-point-fft/scale (twiddle size dst dst-offset src src-offset scale
                            tmp1 tmp2 tmp-offset)
  (declare (ignore twiddle size tmp1 tmp2 tmp-offset)
           (type double-float scale)
           (type complex-sample-array dst src)
           (type index dst-offset src-offset)
           (optimize speed (safety 0)))
  (setf (aref dst dst-offset) (* scale (aref src src-offset)))
  dst)

(defun power-of-two-p (x)
  (and (typep x 'index)
       (= 1 (logcount x))))
