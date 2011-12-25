(in-package "NAPA-FFT")

(deftype half-index ()
  `(and (integer 1)
        (unsigned-byte #.(truncate (integer-length most-positive-fixnum) 2))))

(deftype index ()
  `(unsigned-byte #. (min (1- (integer-length most-positive-fixnum))
                          (integer-length (1- array-dimension-limit)))))

(deftype complex-sample ()
  `(complex double-float))

(deftype complex-sample-array (&optional size)
  `(simple-array complex-sample (,size)))

;; twiddle-factors size
;; dst dst-offset
;; src src-offset
;; [scale]
;; tmp1 tmp2 tmp-offset
(deftype fft-function (&key (scale nil))
  `(function (complex-sample-array index
              complex-sample-array index
              complex-sample-array index
              ,@(and scale '(double-float))
              complex-sample-array complex-sample-array index)
    (values complex-sample-array &optional)))
