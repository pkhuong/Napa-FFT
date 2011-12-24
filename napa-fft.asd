(asdf:defsystem "napa-fft"
  :version "0.0.1"
  :licence "BSD"
  :description "Generates code to execute Fast Fourier Transforms"
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
