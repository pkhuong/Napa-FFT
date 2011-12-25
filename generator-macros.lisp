(in-package "NAPA-FFT")

(defmacro with-fft ((name size direction) &body body)
  (let (generated-name
        names)
    (multiple-value-bind (bodies declares)
        (with-function-bodies
          (setf generated-name (generate-fft-body size direction nil)))
      (setf names (mapcar #'first bodies))
      `(labels (,@bodies
                (,name (twiddle dst dst-offset src src-offset tmp1 tmp2 tmp-offset)
                  (,generated-name twiddle dst dst-offset src src-offset tmp1 tmp2 tmp-offset)))
         (declare ,@declares
                  (inline ,name))
         (flet ,bodies
           (declare ,@declares
                    (ignorable ,@(mapcar (lambda (name) `#',name) names)))
           ,@body)))))

(defmacro with-fft-kernels ((max-fft-size max-transpose-size &key fwd bwd scale)
                            &body body)
  (multiple-value-bind (bodies declares)
      (with-function-bodies
        (loop for i from 1 upto max-fft-size
              for size = (ash 1 i)
              do (when fwd
                   (generate-fft-body size  1 nil)
                   (when scale
                     (generate-fft-body size 1 t)))
                 (when bwd
                   (generate-fft-body size -1 nil)
                   (when scale (generate-fft-body size -1 t))))
        (generate-transpose-body 2 2 nil)
        (generate-transpose-body 2 2 t)
        (generate-inner-transpose-body 2 nil)
        (generate-inner-transpose-body 2 t)
        (loop for i from 2 upto max-transpose-size
              for j = (1- i)
              for size1 = (ash 1 i)
              for size2 = (ash 1 j)
              do (generate-inner-transpose-body size1 nil)
                 (generate-inner-transpose-body size1 t)
                 (generate-transpose-body size1 size1 nil)
                 (generate-transpose-body size1 size1 t)
                 (generate-transpose-body size1 size2 nil)
                 (generate-transpose-body size1 size2 t)
                 (generate-transpose-body size2 size1 nil)
                 (generate-transpose-body size2 size1 t)))
    (let ((ignore `(declare (ignorable ,@(mapcar (lambda (body)
                                                   `#',(first body))
                                                 bodies)))))
      `(labels ,bodies
         (declare ,@declares)
         ,ignore
         (flet ,bodies
           ,ignore
           ,@body)))))

(defun build-blit-routines (fwd-funs bwd-funs fwd/s-funs bwd/s-funs)
  `(progn
     ,@(loop for i upfrom 1
             for fwd in (mapcar (lambda (name)
                                  `(named-lambda ,name
                                       (twiddle size
                                        dst dst-offset
                                        src src-offset
                                        tmp1 tmp2 tmp-offset)
                                     (declare (ignore size))
                                     (,name twiddle dst dst-offset
                                            src src-offset
                                            tmp1 tmp2 tmp-offset)))
                                fwd-funs)
             for bwd in (mapcar (lambda (name)
                                  `(named-lambda ,name
                                       (twiddle size
                                        dst dst-offset
                                        src src-offset
                                        tmp1 tmp2 tmp-offset)
                                     (declare (ignore size))
                                     (,name twiddle dst dst-offset
                                            src src-offset
                                            tmp1 tmp2 tmp-offset)))
                                bwd-funs)
             for fwd/s in (mapcar (lambda (name)
                                    `(named-lambda ,name
                                         (twiddle size
                                          dst dst-offset
                                          src src-offset
                                          scale
                                          tmp1 tmp2 tmp-offset)
                                       (declare (ignore size))
                                       (,name twiddle
                                              dst dst-offset
                                              src src-offset
                                              scale
                                              tmp1 tmp2 tmp-offset)))
                                  fwd/s-funs)
             for bwd/s in (mapcar (lambda (name)
                                    `(named-lambda ,name
                                         (twiddle size
                                          dst dst-offset
                                          src src-offset
                                          scale
                                          tmp1 tmp2 tmp-offset)
                                       (declare (ignore size))
                                       (,name twiddle
                                              dst dst-offset
                                              src src-offset
                                              scale
                                              tmp1 tmp2 tmp-offset)))
                                  bwd/s-funs)
             collect `(setf (aref fwd-fft       ,i) ,fwd
                            (aref bwd-fft       ,i) ,bwd
                            (aref fwd-fft/scale ,i) ,fwd/s
                            (aref bwd-fft/scale ,i) ,bwd/s))))

(defun build-fft-fun (name fft-name vector-name lb-specialized-size
                      &optional (scale nil) scale-vector-name)
  (assert (or (and scale scale-vector-name)
              (not (or scale scale-vector-name))))
  `(,name (twiddle size dst dst-offset src src-offset ,@(and scale '(scale)) tmp1 tmp2 tmp-offset)
            (declare (type complex-sample-array twiddle dst src tmp1 tmp2)
                     ,@(and scale '((type double-float scale)))
                     (type index dst-offset src-offset tmp-offset)
                     (type (and index (integer 1)) size))
            (let ((len (integer-length (1- size))))
              (when (<= len ,lb-specialized-size)
                (return-from ,name (funcall (aref ,(or scale-vector-name vector-name) len)
                                            twiddle size
                                            dst dst-offset
                                            src src-offset
                                            ,@(and scale '(scale))
                                            tmp1 tmp2 tmp-offset)))
              (let ((size1 (ash 1 (truncate len 2)))
                    (size2 (ash 1 (truncate (1+ len) 2))))
                (tran dst dst-offset src src-offset size2 size1)
                (if (<= size1 ,(ash 1 lb-specialized-size))
                    (let ((fun (aref ,vector-name (truncate len 2))))
                      #-xecto-parallel
                      (loop for i of-type index below size by size1 do
                        (let ((tmp-offset (+ tmp-offset i)))
                          (funcall fun
                                   twiddle size1
                                   tmp1 tmp-offset
                                   dst (+ dst-offset i)
                                   tmp2 tmp1 tmp-offset)))
                      #+xecto-parallel
                      (parallel:dotimes (i size2)
                        :wait
                        (declare (type half-index i))
                        (let* ((i          (* i size1))
                               (tmp-offset (+ tmp-offset i)))
                          (funcall fun
                                   twiddle size1
                                   tmp1 tmp-offset
                                   dst (+ dst-offset i)
                                   tmp2 tmp1 tmp-offset))))
                    #-xecto-parallel
                    (loop for i of-type index below size by size1 do
                      (let ((tmp-offset (+ tmp-offset i)))
                        (,fft-name twiddle size1 tmp1 tmp-offset dst (+ dst-offset i) tmp2 tmp1 tmp-offset)))
                    #+xecto-parallel
                    (parallel:dotimes (i size2)
                      :wait
                      (declare (type half-index i))
                      (let* ((i          (* i size1))
                             (tmp-offset (+ tmp-offset i)))
                        (,fft-name twiddle size1 tmp1 tmp-offset dst (+ dst-offset i) tmp2 tmp1 tmp-offset))))
                (tran/twiddle dst dst-offset tmp1 tmp-offset size1 size2 twiddle)
                (if (<= size2 ,(ash 1 lb-specialized-size))
                    (let ((fun (aref ,(or scale-vector-name vector-name)
                                     (truncate (1+ len) 2))))
                      #-xecto-parallel
                      (loop for i of-type index below size by size2 do
                        (let ((tmp-offset (+ tmp-offset i)))
                          (funcall fun
                                   twiddle size2
                                   tmp1 tmp-offset
                                   dst (+ dst-offset i)
                                   ,@(and scale '(scale))
                                   tmp2 tmp1 tmp-offset)))
                      #+xecto-parallel
                      (parallel:dotimes (i size1)
                        :wait
                        (declare (type half-index i))
                        (let* ((i          (* i size2))
                               (tmp-offset (+ tmp-offset i)))
                          (funcall fun
                                   twiddle size2
                                   tmp1 tmp-offset
                                   dst (+ dst-offset i)
                                   ,@(and scale '(scale))
                                   tmp2 tmp1 tmp-offset))))
                    #-xecto-parallel
                    (loop for i of-type index below size by size2 do
                      (let ((tmp-offset (+ tmp-offset i)))
                        (,name twiddle size2
                               tmp1 tmp-offset
                               dst (+ dst-offset i)
                               ,@(and scale '(scale))
                               tmp2 tmp1 tmp-offset)))
                    #+xecto-parallel
                    (parallel:dotimes (i size1)
                        :wait
                      (declare (type half-index i))
                      (let* ((i          (* i size2))
                             (tmp-offset (+ tmp-offset i)))
                        (,name twiddle size2
                               tmp1 tmp-offset
                               dst (+ dst-offset i)
                               ,@(and scale '(scale))
                               tmp2 tmp1 tmp-offset))))
                (tran dst dst-offset tmp1 tmp-offset size2 size1)))))

(defun build-transpose-funs (outer-name inner-name twiddle half-size base-case)
  `((,inner-name (dst dst-offset dst-stride src src-offset src-stride size
                      ,@(and twiddle '(twiddle twiddle-offset)))
      (declare (type complex-sample-array dst src ,@(and twiddle '(twiddle)))
               (type index dst-offset src-offset ,@(and twiddle '(twiddle-offset)))
               (type half-index dst-stride src-stride size))
      (cond ((= size ,half-size)
             ;; other cases are picked off by specialised FFTs
             (,base-case
              dst dst-offset dst-stride
              src src-offset src-stride
              ,@(and twiddle '(twiddle twiddle-offset))))
            (t
             (let* ((size/2          (truncate size 2))
                    (long-dst-stride (* size/2 dst-stride))
                    (long-src-stride (* size/2 src-stride)))
               (declare (type index long-dst-stride long-src-stride))
               #-xecto-parallel
               (unrolled-for ((i 4))
                 (let* ((short (- (logand i 1)))
                        (long  (- (logand i 2)))
                        (dst-delta (+ (logand size/2 short)
                                      (logand long-dst-stride long)))
                        (src-delta (+ (logand long-src-stride short)
                                      (logand size/2 long))))
                   (declare (type index dst-delta src-delta))
                   (,inner-name dst (+ dst-offset dst-delta)
                                dst-stride
                                src (+ src-offset src-delta)
                                src-stride
                                size/2
                                ,@(and twiddle '(twiddle (+ twiddle-offset dst-delta))))))
               #+xecto-parallel
               (parallel:let ((aa (,inner-name dst dst-offset
                                               dst-stride
                                               src src-offset
                                               src-stride
                                               size/2
                                               ,@(and twiddle '(twiddle twiddle-offset))))
                              (ab (,inner-name dst (+ dst-offset size/2)
                                               dst-stride
                                               src (+ src-offset long-src-stride)
                                               src-stride
                                               size/2
                                               ,@(and twiddle '(twiddle (+ twiddle-offset size/2)))))
                              (ba (,inner-name dst (+ dst-offset long-dst-stride)
                                               dst-stride
                                               src (+ src-offset size/2)
                                               src-stride
                                               size/2
                                               ,@(and twiddle '(twiddle (+ twiddle-offset long-dst-stride)))))
                              (bb (,inner-name dst (+ dst-offset long-dst-stride size/2)
                                               dst-stride
                                               src (+ src-offset long-src-stride size/2)
                                               src-stride
                                               size/2
                                               ,@(and twiddle '(twiddle (+ twiddle-offset long-dst-stride size/2)))))
                              (:parallel (>= (* size/2 size/2) 1024)))
                 (declare (ignore aa ab ba bb)))))))
    (,outer-name (dst dst-offset src src-offset size1 size2
                      ,@(and twiddle '(twiddle &aux (twiddle-base (+ (* size1 size2)
                                                                     +factor-bias+)))))
      (declare (type complex-sample-array dst src ,@(and twiddle '(twiddle)))
               (type index dst-offset src-offset ,@(and twiddle '(twiddle-base)))
               (type half-index size1 size2))
      (cond ((= size1 size2)
             (,inner-name dst dst-offset size2 src src-offset size1 size1
                          ,@(and twiddle '(twiddle twiddle-base))))
            ((< size1 size2)
             (let* ((size  size1)
                    (block (* size size)))
               (#-xecto-parallel let
                #+xecto-parallel parallel:let
                ((a (,inner-name dst dst-offset size2
                                 src src-offset size1
                                 size
                                 ,@(and twiddle '(twiddle twiddle-base))))
                 (b (,inner-name dst (+ size  dst-offset) size2
                                 src (+ block src-offset) size1
                                 size
                                 ,@(and twiddle '(twiddle (+ size twiddle-base)))))
                 #+xecto-parallel (:parallel (>= (* size size) 1024)))
                (declare (ignore a b)))))
            (t
             (let* ((size  size2)
                    (block (* size size)))
               (#-xecto-parallel let
                #+xecto-parallel parallel:let
                ((a (,inner-name dst dst-offset size2
                                 src src-offset size1
                                 size
                                 ,@(and twiddle '(twiddle twiddle-base))))
                 (b (,inner-name dst (+ block dst-offset) size2
                                 src (+ size  src-offset) size1
                                 size
                                 ,@(and twiddle '(twiddle (+ block twiddle-base)))))
                 #+xecto-parallel (:parallel (>= (* size size) 1024)))
                 (declare (ignore a b)))))))))

(defmacro build-fft-routine-vectors (lb-specialized-size lb-max-size)
  (when (oddp lb-specialized-size)
    (incf lb-specialized-size))
  (setf lb-max-size (max lb-specialized-size lb-max-size))
  (let* ((count (1+ lb-max-size))
         (fwd-funs   (loop for i from 1 upto lb-specialized-size
                           collect (symbolicate 'fft (ash 1 i) :fwd)))
         (bwd-funs   (loop for i from 1 upto lb-specialized-size
                           collect (symbolicate 'fft (ash 1 i) :bwd)))
         (fwd/s-funs (loop for i from 1 upto lb-specialized-size
                           collect (symbolicate 'fft (ash 1 i) :fwd :scale)))
         (bwd/s-funs (loop for i from 1 upto lb-specialized-size
                           collect (symbolicate 'fft (ash 1 i) :bwd :scale)))
         (half-size  (ash 1 (truncate lb-specialized-size 2)))
         (%tran-base-case (symbolicate '%transpose half-size))
         (%tran/twiddle-base-case (symbolicate '%transpose half-size :twiddle)))
    `(let ((fwd-fft         (make-array ,count))
           (bwd-fft         (make-array ,count))
           (fwd-fft/scale   (make-array ,count))
           (bwd-fft/scale   (make-array ,count)))
       (declare (type (simple-array fft-function (,count))
                      fwd-fft bwd-fft)
                (type (simple-array (fft-function :scale t) (,count))
                      fwd-fft/scale bwd-fft/scale)
                (optimize speed (safety 0)))
       (with-fft-kernels (,lb-specialized-size ,(truncate lb-specialized-size 2)
                          :fwd t :bwd t :scale t)
         (declare (inline ,@fwd-funs ,@bwd-funs
                          ,@fwd/s-funs ,@bwd/s-funs
                          ,%tran-base-case ,%tran/twiddle-base-case))
         (setf (aref fwd-fft 0) #'one-point-fft)
         (setf (aref bwd-fft 0) #'one-point-fft)
         (setf (aref fwd-fft/scale 0) #'one-point-fft/scale)
         (setf (aref bwd-fft/scale 0) #'one-point-fft/scale)
         ,(build-blit-routines fwd-funs bwd-funs fwd/s-funs bwd/s-funs)
         (labels (,(build-fft-fun 'fft/fwd       'fft/fwd 'fwd-fft lb-specialized-size)
                  ,(build-fft-fun 'fft/fwd-scale 'fft/fwd 'fwd-fft lb-specialized-size
                     t 'fwd-fft/scale)
                  ,(build-fft-fun 'fft/bwd       'fft/bwd 'bwd-fft lb-specialized-size)
                  ,(build-fft-fun 'fft/bwd-scale 'fft/bwd 'bwd-fft lb-specialized-size
                     t 'bwd-fft/scale)
                  ,@(build-transpose-funs 'tran '%tran
                      nil half-size
                      %tran-base-case)
                  ,@(build-transpose-funs 'tran/twiddle '%tran/twiddle-base-case
                      t half-size
                      %tran/twiddle-base-case))
           (declare (inline tran tran/twiddle))
           ,@(loop for vec in '(fwd-fft fwd-fft/scale bwd-fft bwd-fft/scale)
                   for fun in '(fft/fwd fft/fwd-scale fft/bwd fft/bwd-scale)
                   collect `(fill ,vec #',fun :start ,(1+ lb-specialized-size)))
           (values fwd-fft
                   fwd-fft/scale
                   bwd-fft
                   bwd-fft/scale
                   #'fft/fwd
                   #'fft/fwd-scale
                   #'fft/bwd
                   #'fft/bwd-scale))))))
