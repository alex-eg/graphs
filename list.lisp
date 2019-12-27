(defpackage :cl-math.list
  (:use :cl)
  (:export
   :find-all
   :remove-all))

(in-package :cl-math.list)

(defun find-all (item list &key (key #'identity) (test #'=))
  "Collect all elements conforming to `test' in a list"
  (let (res)
    (dolist (e list res)
      (if (funcall test (funcall key e) item)
        (push e res)))))

(defun remove-all (item-list list)
  "Remove each element of `item-list' from `list'"
  (reduce (lambda (a e) (remove e a))
          item-list :initial-value list))
