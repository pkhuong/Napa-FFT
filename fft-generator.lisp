(in-package "NAPA-FFT")

(defvar *fft-inline* 4)

(defun generate-fft/2 (scale &optional (steps 1) (stepd 1))
  `(lambda (twiddle dst dst-offset src src-offset ,@(and scale '(scale)))
     (declare (type index dst-offset src-offset)
              (type complex-sample-array dst src)
              ,@(and scale '((type double-float scale)))
              (ignore twiddle))
     (let ((s0 (aref src     src-offset))
           (s1 (aref src (+ src-offset ,steps))))
       ,(if scale
            `(setf (aref dst     dst-offset)        (* scale (+ s0 s1))
                   (aref dst (+ dst-offset ,stepd)) (* scale (- s0 s1)))
            `(setf (aref dst     dst-offset)        (+ s0 s1)
                   (aref dst (+ dst-offset ,stepd)) (- s0 s1)))
       dst)))

(defun generate-fft/4 (direction scale &optional (steps 1) (stepd 1))
  (flet ((scalify (form)
           (if scale
               `(* scale ,form)
               form)))
    `(lambda (twiddle dst dst-offset src src-offset ,@(and scale '(scale)))
       (declare (type index dst-offset src-offset)
                (type complex-sample-array dst src)
                ,@(and scale '((type double-float scale)))
                (ignore twiddle))
       ,(if (plusp direction)
            `(let* ((s0 (aref src                 src-offset))
                    (s2 (aref src (+ src-offset ,(* 2 steps))))
                    (s0+s2 (+ s0 s2))
                    (s0-s2 (- s0 s2))
                    (s1 (aref src (+      src-offset  ,steps)))
                    (s3 (aref src (+ src-offset ,(* 3 steps))))
                    (s1+s3 (+ s1 s3))
                    (s1-s3 (- (mul+i s1)
                              (mul+i s3))))
               (setf (aref dst                 dst-offset)  ,(scalify `(+ s0+s2 s1+s3))
                     (aref dst (+      dst-offset  ,stepd)) ,(scalify `(- s0-s2 s1-s3))
                     (aref dst (+ dst-offset ,(* 2 stepd))) ,(scalify `(- s0+s2 s1+s3))
                     (aref dst (+ dst-offset ,(* 3 stepd))) ,(scalify `(+ s0-s2 s1-s3)))
               dst
               #+nil
               (setf (aref dst      dst-offset)  (+ s0 s1 s2 s3)
                     (aref dst (+ 1 dst-offset)) (+ s0
                                                    (mul-i s1)
                                                    (- s2)
                                                    (mul+i s3))
                     (aref dst (+ 2 dst-offset)) (+ s0
                                                    (- s1)
                                                    (+ s2)
                                                    (- s3))
                     (aref dst (+ 3 dst-offset)) (+ s0
                                                    (mul+i s1)
                                                    (- s2)
                                                    (mul-i s3))))
            `(let* ((s0 (aref src      src-offset))
                    (s2 (aref src (+ src-offset ,(* 2 steps))))
                    (s0+s2 (+ s0 s2))
                    (s0-s2 (- s0 s2))
                    (s1 (aref src (+      src-offset  ,steps)))
                    (s3 (aref src (+ src-offset ,(* 3 steps))))
                    (s1+s3 (+ s1 s3))
                    (s1-s3 (- (mul+i s1)
                              (mul+i s3))))
               (setf (aref dst                 dst-offset)  ,(scalify `(+ s0+s2 s1+s3))
                     (aref dst (+      dst-offset  ,stepd)) ,(scalify `(+ s0-s2 s1-s3))
                     (aref dst (+ dst-offset ,(* 2 stepd))) ,(scalify `(- s0+s2 s1+s3))
                     (aref dst (+ dst-offset ,(* 3 stepd))) ,(scalify `(- s0-s2 s1-s3)))
               dst)))))

(defun generate-fft/8 (direction scale)
  ;; size1: 2
  ;; size2: 4
  (let ((fft-size1 (generate-fft/2 nil 4))
        (fft-size2 (generate-fft/4 direction scale 1 2)))
    `(lambda (twiddle dst dst-offset src src-offset ,@(and scale '(scale)) tmp tmp-offset)
       (declare (type complex-sample-array dst src twiddle)
                ,@(and scale '((type double-float scale))))
       ,@(loop for i below 4
               collect `(,fft-size1 twiddle
                                    dst (+ dst-offset ,(* 2 i))
                                    src (+ src-offset ,i)))
       (,(generate-transpose-body 2 4 direction)
        tmp tmp-offset dst dst-offset twiddle)
       (,@(if scale
              `(let ((scale2 0d0))
                 (declare (type double-float scale2))
                 (setf scale2 scale))
              '(progn))
        ,@(loop for i below 2
                collect `(,fft-size2 twiddle
                                     dst (+ dst-offset ,i)
                                     tmp (+ tmp-offset ,(* 4 i))
                                     ,@(and scale '(scale2)))))
       dst)))

(defun generate-fft/16 (direction scale)
  (let ((fft-size1 (generate-fft/4 direction nil 4))
        (fft-size2 (generate-fft/4 direction scale 1 4)))
    `(lambda (twiddle dst dst-offset src src-offset
              ,@(and scale '(scale)) tmp tmp-offset)
       (declare (type complex-sample-array twiddle dst src tmp)
                (type index dst-offset src-offset tmp-offset)
                ,@(and scale '((type double-float scale))))
       (prog ((i  3)
              (j 12))
          (declare (type index i j))
          (go body)
        loop
          (decf i)
          (decf j 4)
        body
          (,fft-size1 twiddle
                      dst (+ dst-offset j)
                      src (+ src-offset i))
          (when (zerop j)
            (return))
          (go loop))
       (,(generate-transpose-body 4 4 direction)
        tmp tmp-offset dst dst-offset twiddle)
       (prog ((i  3)
              (j 12)
              ,@(and scale '((scale2 0d0))))
          (declare (type index i j)
                   ,@(and scale '((type double-float scale2))))
          ,@(and scale '((setf scale2 scale)))
          (go body)
        loop
          (decf i)
          (decf j 4)
        body
          (,fft-size2 twiddle
                      dst (+ dst-offset i)
                      tmp (+ tmp-offset j)
                      ,@(and scale '(scale2)))
          (when (zerop i)
            (return))
          (go loop))
       dst)))

(defun generate-fft-body (size direction scale)
  (flet ((gen ()
           (values
            `((declare (ignorable tmp1 tmp2 tmp-offset)
                       (type complex-sample-array twiddle dst src tmp1 tmp2)
                       (type index dst-offset src-offset tmp-offset)
                       (ignorable twiddle)
                       ,@(and scale `((type double-float scale))))
              ,(case size
                 (2
                  `(,(generate-fft/2 scale) twiddle dst dst-offset src src-offset
                    ,@(and scale '(scale))))
                 (4
                  `(,(generate-fft/4 direction scale) twiddle dst dst-offset src src-offset
                    ,@(and scale '(scale))))
                 (8
                  `(,(generate-fft/8 direction scale) twiddle dst dst-offset src src-offset
                    ,@(and scale '(scale)) tmp1 tmp-offset))
                 (16
                  `(,(generate-fft/16 direction scale) twiddle dst dst-offset src src-offset
                    ,@(and scale '(scale)) tmp1 tmp-offset))
                 (t
                  (let* ((lb    (integer-length (1- size)))
                         (size1 (ash 1 (truncate lb 2)))
                         (size2 (the integer (/ size size1)))
                         (unroll-1 (and (> size1 *fft-inline*)
                                        (<= size2 4)))
                         (unroll-2 (and (> size2 *fft-inline*)
                                        (<= size1 4))))
                    `(progn
                       (,(generate-transpose-body size2 size1) dst dst-offset src src-offset)
                       ,@(if unroll-1
                             (loop with op = (generate-fft-body size1 direction nil)
                                   for i from 0 below size by size1
                                   collect `(,op twiddle
                                                 tmp1 (+ tmp-offset ,i)
                                                 dst  (+ dst-offset ,i)
                                                 tmp2 tmp1 (+ tmp-offset ,i)))
                             `((loop for i of-type index from 0 below ,size by ,size1
                                     do (let ((tmp-offset (+ tmp-offset i)))
                                          (,(generate-fft-body size1 direction nil)
                                           twiddle
                                           tmp1 tmp-offset
                                           dst  (+ dst-offset i)
                                           tmp2 tmp1 tmp-offset)))))
                       (,(generate-transpose-body size1 size2 direction)
                        dst dst-offset tmp1 tmp-offset twiddle)
                       (,@(if scale
                              '(let ((scale2 0d0))
                                (declare (type double-float scale2))
                                (setf scale2 scale))
                              '(progn))
                        ,@(if unroll-2
                              (loop with op = (generate-fft-body size2 direction scale)
                                    for i from 0 below size by size2
                                    collect `(,op twiddle
                                                  tmp1 (+ tmp-offset ,i)
                                                  dst (+ dst-offset ,i)
                                                  ,@(and scale '(scale2))
                                                  tmp2 tmp1 (+ tmp-offset ,i)))
                              `((loop for i of-type index from 0 below ,size by ,size2
                                      do (let ((tmp-offset (+ tmp-offset i)))
                                           (,(generate-fft-body size2 direction scale)
                                            twiddle
                                            tmp1 tmp-offset
                                            dst (+ dst-offset i)
                                            ,@(and scale '(scale2))
                                            tmp2 tmp1 tmp-offset))))))
                       (,(generate-transpose-body size2 size1) dst dst-offset tmp1 tmp-offset)
                       dst)))))
            (and (<= size *fft-inline*)
                 `((inline .self.))))))
    (if scale
        (ensure-body (fft size (if (plusp direction) :fwd :bwd) :scale)
            (twiddle dst dst-offset src src-offset scale tmp1 tmp2 tmp-offset)
          (gen))
        (ensure-body (fft size (if (plusp direction) :fwd :bwd))
            (twiddle dst dst-offset src src-offset tmp1 tmp2 tmp-offset)
          (gen)))))
