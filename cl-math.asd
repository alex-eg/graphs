(defpackage :cl-math-asd
  (:use :cl :asdf))

(in-package :cl-math-asd)

(asdf:defsystem cl-math
  :components
  ((:file "io")
   (:file "set")
   (:file "matrix")
   (:file "polynomial")
   (:file "graph")))
