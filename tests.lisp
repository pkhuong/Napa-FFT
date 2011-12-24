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
         (repeats  128)
         (inner-loop (if (> size (ash 1 16)) 1 10))
         (normal-times (make-array repeats))
         (fft-times    (make-array repeats))
         (factors  *fwd-factors*)
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

(defun test-fft-vector (vector &optional (max (length vector)))
  (loop for i from 1 below max
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

