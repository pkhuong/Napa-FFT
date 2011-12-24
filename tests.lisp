(in-package "NAPA-FFT")

(defun quantile (vector quantile)
  (let ((vector (sort vector #'<)))
    (aref vector (truncate (* (length vector) quantile)))))

(defun test-ffts (fft size)
  (let* ((source (make-array size :element-type 'complex-sample :initial-element (complex 1d0)))
         (dest1  (make-array size :element-type 'complex-sample :initial-element (complex 1d0)))
         (dest2  (make-array size :element-type 'complex-sample :initial-element (complex 1d0)))
         (tmp1   (make-array size :element-type 'complex-sample :initial-element (complex 1d0)))
         (tmp2   (make-array size :element-type 'complex-sample :initial-element (complex 1d0)))
         (instance (bordeaux-fft::make-fourier-instance size 1))
         (repeats  10)
         (inner-loop (if (> size (ash 1 16)) 1 10))
         (normal-times (make-array repeats))
         (fft-times    (make-array repeats))
         (factors    (ensure-factor-size *fwd-factors* size 1))
         (scale   (* 1d0 inner-loop (* size (integer-length (1- size))))))
    (dotimes (i repeats (values (/ (quantile normal-times 1d-1)
                                   scale)
                                (/ (quantile fft-times 1d-1)
                                   scale)))
      (dotimes (j size)
        (setf (aref source j) (- (complex (random 2d0) (random 2d0))
                                 #c(1d0 1d0))))
      (setf (aref normal-times i)
            (nth-value
             1 (sb-vm::with-cycle-counter
                 (dotimes (i inner-loop)
                   (bordeaux-fft::fft-common instance source dest1)))))
      (setf (aref fft-times i)
            (nth-value
             1 (sb-vm::with-cycle-counter
                 (dotimes (i inner-loop)
                   (funcall fft factors size dest2 0 source 0 tmp1 tmp2 0)))))
      #+nil
      (map-into dest2 (lambda (x) (* 2 x)) dest2)
      (unless (every (lambda (x y)
                       (< (abs (- x y)) 1d-3))
                     dest1 dest2)
        (break "Discrepancy at ~A: ~A ~A~%" i dest1 dest2)))))

(defun test-fft-vector (vector &optional (max (length vector)) (min 1))
  (loop for i from min below max
        for len = (ash 1 i)
        do (let ((fun (aref vector i)))
             (sb-ext:gc :full t)
             (format t "size: ~A~%" i)
             (format t "   ~{~A ~}~%" (multiple-value-bind (slow fast)
                                          (test-ffts fun len)
                                        (list (/ slow fast) fast slow))))))

#+nil
(defun test-all-sizes (max)
  (loop for i from 1 upto max
        for len = (ash 1 i)
        do (let ((fun (let ((*error-output* (make-broadcast-stream)))
                        (compile nil `(lambda (dst src tmp1 tmp2)
                                        (declare (complex-sample-array dst src tmp1 tmp2)
                                                 (optimize speed (safety 0)))
                                        (recursive-fft dst src ,len tmp1 tmp2)
                                        #+nil
                                        (with-fft (fft ,len 1)
                                          (fft dst 0 src 0 tmp1 tmp2 0)))))))
             (format t "size: ~A~%" len)
             (format t "   ~{~A ~}~%" (multiple-value-bind (slow fast)
                                          (test-ffts fun len)
                                        (list (/ slow fast) fast slow))))))

#||
Parallel:
* (parallel-future:with-context (8) (test-fft-vector *forward* 27 20))
size: 20
   4.753233346171666d0 2.8306907653808593d0 13.454933738708496d0
size: 21
   4.915500455687723d0 2.792923791067941d0 13.728618167695545d0
size: 22
   5.336004697772222d0 2.6783089421012183d0 14.291469097137451d0
size: 23
   5.890097446393686d0 2.591226691785066d0 15.262577720310377d0
size: 24
   6.262575770053327d0 2.4720080445210137d0 15.481137682994207d0
size: 25
   6.628124577258128d0 2.3695152401924133d0 15.705442199707031d0
size: 26
   7.147543084293575d0 2.246080761918655d0 16.053959016616528d0 
NIL


||#
