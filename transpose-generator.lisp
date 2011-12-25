(in-package "NAPA-FFT")

(defvar *transpose-unroll* 64)

(defun generate-unrolled-transpose-body (size twiddle)
  `((unrolled-for ((i ,size)
                   (src-offset :offset src-offset
                               :type index)
                   (dst-offset :offset dst-offset :stride dst-stride
                               :type index)
                   ,@(and twiddle `((twiddle-offset :offset twiddle-offset
                                                    :stride dst-stride
                                                    :type index))))
      (unrolled-for ((j ,size)
                     (src-offset :offset src-offset :stride src-stride
                                 :type index)
                     (dst-offset :offset dst-offset
                                 :type index)
                     ,@(and twiddle `((twiddle-offset :offset twiddle-offset
                                                      :type index))))
        (setf (aref dst dst-offset) ,(if twiddle
                                         `(* (aref twiddle twiddle-offset)
                                             (aref src src-offset))
                                         `(aref src src-offset)))))
    dst))

(defun generate-inner-transpose-body (size &optional twiddle)
  (flet ((generate-body ()
           `((declare (type complex-sample-array dst src ,@(and twiddle '(twiddle)))
                      (type index dst-offset src-offset  ,@(and twiddle '(twiddle-offset)))
                      (type half-index dst-stride src-stride))
             ,@(if (<= (* size size) *transpose-unroll*)
                   (generate-unrolled-transpose-body size twiddle)
                   (let* ((size/2 (truncate size 2))
                          (op     (generate-inner-transpose-body size/2 twiddle)))
                     `((let ((long-dst-stride (* dst-stride ,size/2))
                             (long-src-stride (* src-stride ,size/2)))
                         (declare (type index long-dst-stride long-src-stride))
                         (,@(if  (<= (* size/2 size/2) *transpose-unroll*)
                                 `(loop for i below 4 do)
                                 `(unrolled-for ((i 4))))
                          (let* ((short (- (logand i 1)))
                                 (long  (- (logand i 2)))
                                 (dst-delta (+ (logand ,size/2 short)
                                               (logand long-dst-stride long)))
                                 (src-delta (+ (logand long-src-stride short)
                                               (logand ,size/2 long))))
                            (declare (type index dst-delta src-delta))
                            ;; masks are ok: size/2 >= 2
                            (,op dst (+ dst-offset dst-delta) dst-stride
                                 src (+ src-offset src-delta) src-stride
                                 ,@(and twiddle
                                        `(twiddle (+ twiddle-offset dst-delta))))))
                         dst)))))))
    (if (not twiddle)
        (ensure-body (%transpose size)
            (dst dst-offset dst-stride src src-offset src-stride)
          (generate-body))
        (ensure-body (%transpose size :twiddle)
            (dst dst-offset dst-stride src src-offset src-stride
             twiddle twiddle-offset)
          (generate-body)))))

(defun generate-transpose-body (size1 size2 &optional twiddle)
  (declare (type (member nil -1 1 t) twiddle))
  (flet ((gen (&aux (total-size (* size1 size2))
                    (twiddle-base (+ +factor-bias+ total-size)))
           (values
            `((declare (type complex-sample-array dst src ,@(and (eql twiddle t) '(twiddle)))
                       ,@(and twiddle '((ignorable twiddle)))
                       (type index dst-offset src-offset))
              ,@(cond ((and (<= (* size1 size2) *transpose-unroll*)
                            (member twiddle '(nil -1 1)))
                       `((setf
                          ,@(loop for i below size1
                                  nconc (loop
                                          for j below size2
                                          for dst-index = (+ (* i size2) j)
                                          for src = `(aref src (+ ,(+ (* j size1) i) src-offset))
                                          collect `(aref dst (+ ,dst-index dst-offset))
                                          collect (if twiddle
                                                      (let ((factor (/ (mod (* i j -1 twiddle)
                                                                            total-size)
                                                                       total-size))
                                                            (coef   `(aref twiddle ,(+ twiddle-base
                                                                                       dst-index))))
                                                        (case factor
                                                          (0 src)
                                                          (1/2
                                                           `(- ,src))
                                                          (1/4
                                                           `(mul+i ,src))
                                                          (3/4
                                                           `(mul-i ,src))
                                                          ((1/8 5/8)
                                                           `(mul+/-sqrt+i ,src (realpart ,coef)))
                                                          ((3/8 7/8)
                                                           `(mul+/-sqrt-i ,src (realpart ,coef)))
                                                          (t
                                                           `(* ,coef ,src))))
                                                      src))))
                         dst))
                      ((= size1 size2)
                       `((,(generate-inner-transpose-body size1 twiddle)
                           dst dst-offset ,size2
                           src src-offset ,size1
                           ,@(and twiddle `(twiddle ,twiddle-base)))))
                      ((< size1 size2)
                       (assert (= (* 2 size1) size2))
                       (let* ((size  size1)
                              (block (* size size))
                              (op    (generate-inner-transpose-body size twiddle)))
                         (cond ((> block *transpose-unroll*)
                                `((,op dst dst-offset ,size2
                                       src src-offset ,size1
                                       ,@(and twiddle `(twiddle ,twiddle-base)))
                                  (,op dst (+ dst-offset  ,size) ,size2
                                       src (+ src-offset ,block) ,size1
                                       ,@(and twiddle `(twiddle ,(+ twiddle-base size))))))
                               (twiddle
                                `((loop repeat 2
                                        for offset of-type index
                                          from 0          by ,size
                                        for src-offset of-type index
                                          from src-offset by ,block
                                        do (,op dst     (+ offset dst-offset) ,size2
                                                src     src-offset            ,size1
                                                twiddle (+ offset ,twiddle-base)))))
                               (t
                                `((loop repeat 2
                                        for dst-offset of-type index
                                          from dst-offset by ,size
                                        for src-offset of-type index
                                          from src-offset by ,block
                                        do (,op dst dst-offset ,size2
                                                src src-offset ,size1)))))))
                      (t
                       (assert (= (* 2 size2) size1))
                       (let* ((size  size2)
                              (block (* size size))
                              (op    (generate-inner-transpose-body size twiddle)))
                         (cond ((> block *transpose-unroll*)
                                `((,op dst dst-offset ,size2
                                       src src-offset ,size1
                                       ,@(and twiddle `(twiddle ,twiddle-base)))
                                  (,op dst (+ ,block dst-offset) ,size2
                                       src (+ ,size  src-offset) ,size1
                                       ,@(and twiddle `(twiddle ,(+ twiddle-base block))))))
                               (twiddle
                                `((loop repeat 2
                                        for offset of-type index
                                          from 0          by ,block
                                        for src-offset of-type index
                                          from src-offset by ,size
                                        do (,op dst     (+ offset dst-offset) ,size2
                                                src     src-offset            ,size1
                                                twiddle (+ offset ,twiddle-base)))))
                               (t
                                `((loop repeat 2
                                        for dst-offset of-type index
                                          from dst-offset by ,block
                                        for src-offset of-type index
                                          from src-offset by ,size
                                        do (,op dst dst-offset ,size2
                                                src src-offset ,size1)))))))))
            (and (or (<= (* size1 size2) *transpose-unroll*)
                     (= size1 size2))
                 `((inline .self.))))))
    (ecase twiddle
      ((-1 1)
       (ensure-body (transpose size1 size2 (if (plusp twiddle)
                                               :fwd :bwd))
           (dst dst-offset src src-offset twiddle)
          (gen)))
      ((t)
       (ensure-body (transpose size1 size2 :any)
           (dst dst-offset src src-offset twiddle)
         (gen)))
      ((nil)
       (ensure-body (transpose size1 size2) (dst dst-offset src src-offset)
         (gen))))))
