(defpackage :cl-math-asd
  (:use :cl :asdf))

(in-package :cl-math-asd)

(asdf:defsystem cl-math
  :components
  ((:file "set")
   (:file "polynomial")
   (:file "list")
   (:file "matrix")
   (:file "graph")
   (:file "io")))
