(asdf:defsystem "napa-fft"
  :version "0.0.1"
  :licence "BSD"
  :description "Fast Fourier Transforms in Common Lisp based on a code generator"
  :serial t
  :components ((:file "package")
               (:file "types")
               (:file "generator-support")
               (:file "runtime-support")
               (:file "transpose-generator")
               (:file "fft-generator")
               (:file "generator-macros")
               (:file "interface")
               (:file "windowing")))
