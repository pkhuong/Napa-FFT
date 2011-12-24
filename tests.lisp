(in-package "NAPA-FFT")

(defun test-ffts (fft size)
  (let ((source (make-array size :element-type 'complex-sample :initial-element (complex 1d0)))
        (dest1  (make-array size :element-type 'complex-sample :initial-element (complex 1d0)))
        (dest2  (make-array size :element-type 'complex-sample :initial-element (complex 1d0)))
        (tmp1   (make-array size :element-type 'complex-sample :initial-element (complex 1d0)))
        (tmp2   (make-array size :element-type 'complex-sample :initial-element (complex 1d0)))
        (instance (bordeaux-fft::make-fourier-instance size 1))
        (normal-times 0d0)
        (fft-times   0d0)
        (repeats  128))
    (dotimes (i repeats (values (/ normal-times repeats (* size (integer-length (1- size))))
                                (/ fft-times    repeats (* size (integer-length (1- size))))))
      (dotimes (j size)
        (setf (aref source j) (- (complex (random 2d0) (random 2d0))
                                 #c(1d0 1d0))))
      (incf normal-times
            (nth-value
             1 (sb-vm::with-cycle-counter
                 (bordeaux-fft::fft-common instance source dest1))))
      (incf fft-times
            (nth-value
             1 (sb-vm::with-cycle-counter
                 (funcall fft *fwd-factors* size dest2 0 source 0 tmp1 tmp2 0))))
      #+nil
      (map-into dest2 (lambda (x) (* 2 x)) dest2)
      (unless (every (lambda (x y)
                       (< (abs (- x y)) 1d-3))
                     dest1 dest2)
        (break "Discrepancy at ~A: ~A ~A~%" i dest1 dest2)))))

(defun test-fft-vector (vector)
  (loop for i from 1 below (length vector)
        for len = (ash 1 i)
        do (let ((fun (aref vector i)))
             (sb-ext:gc :full t)
             (format t "size: ~A~%" len)
             (format t "   ~{~A ~}~%" (multiple-value-bind (slow fast)
                                          (test-ffts fun len)
                                        (list (/ slow fast) fast slow))))))

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

