(defpackage :cl-math.set
  (:use :cl)
  (:export
   :subset-find
   :subset-remove
   :subset-merge))

(in-package :cl-math.set)

;;; Sets use list as container

(defun subset-find (set v)
  "Finds a subset containing and element v"
  (find v set :test (lambda (v set)
                      (find v set))))

(defun subset-remove (set subset)
  (delete subset set :test #'eq))

(defun subset-merge (set subset-1 subset-2)
  (nconc subset-1 subset-2)
  (subset-remove set subset-2))
