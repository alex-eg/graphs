(defpackage :cl-math.io
  (:use
   :cl
   :cl-math.polynomial)
  (:import-from :cl-math.graph :set-edge)
  (:import-from :cl-math.matrix :make-square-matrix)
  (:export
   ;; Misc
   :println

   ;; Graph I/O
   ;; */contest functions expect graphs in contest format:
   ;; v e
   ;; v1 v2 w1
   ;; ...
   ;; where v is number of vertices, e is number of edges,
   ;; followed by e lines (1 per edge), where v1 and v2
   ;; are vertex numbers and w1 is edge weight
   :read-weighted-graph-from-string/contest
   :read-weighted-graph/contest

   ;; Matrix
   :print-matrix

   ;; Polynomials
   :print-polynomial))

(in-package :cl-math.io)

(defun println (obj &optional name)
  "It's like print, but with newline at the end.
Returns its argument"
  (if name
    (progn (princ name)
           (princ " ")))
  (princ obj)
  (terpri)
  obj)

(defun string-to-list (string)
  "Read from string as reader would do, and return it as list.
Example: (string-to-list \"1 2 3 wut quux\") -> '(1 2 3 wut quux)"
  (with-input-from-string (s string)
    (loop
      :for c := (read s nil)
      :while c :collect c)))

;;; Graph I/O

(defun read-weighted-graph/contest ()
  "Read graph from *standard-input*, return it in matrix form.
Worth noting that contest format is usually 1-based, so we subtract 1 from
indexes"
  (let* ((nodes (read))
         (edges (read))
         (m (make-square-matrix nodes)))
    (loop
      :repeat edges
      :do (apply
           #'set-edge
           (append
            (list m)
            (let ((s (string-to-list (read-line))))
              (setf (car s) (1- (car s)))
              (setf (cadr s) (1- (cadr s)))
              s))))
    m))

(defun read-weighted-graph-from-string/contest (string)
  (with-input-from-string (*standard-input* string)
    (read-weighted-graph/contest)))

(defun print-matrix (m &key (print-function #'princ))
  "Prints matrix and returns it.
Useful for in-place printing at any step, since the argument matrix is
returned"
  (let ((row-length (array-dimension m 1)))
    (loop
      :for i :below (array-total-size m)
      :do (if (zerop (mod i row-length))
            (terpri)
            (princ #\Tab))
          (funcall print-function (row-major-aref m i))))
  (terpri)
  m)

;;; Polynomials

(defun number-superscript (number)
  (labels ((reverse-superscript (number)
             (let ((super-scripts (list
                                   #\SUPERSCRIPT_ZERO
                                   #\SUPERSCRIPT_ONE
                                   #\SUPERSCRIPT_TWO
                                   #\SUPERSCRIPT_THREE
                                   #\SUPERSCRIPT_FOUR
                                   #\SUPERSCRIPT_FIVE
                                   #\SUPERSCRIPT_SIX
                                   #\SUPERSCRIPT_SEVEN
                                   #\SUPERSCRIPT_EIGHT
                                   #\SUPERSCRIPT_NINE)))
               (if (= 0 number) nil
                 (cons (nth (mod number 10) super-scripts)
                       (reverse-superscript (floor number 10)))))))
    (coerce (nreverse (reverse-superscript number)) 'string)))

(defun print-monomial (m &key print-sign)
  (let ((pwr (pwr m))
        (cff (cff m)))
    (cond ((= cff 0) (princ 0) (return-from print-monomial))
          ((and (= cff 1) print-sign) (princ "+"))
          ((= cff -1) (princ "-"))
          (t (if (and (> cff 1) print-sign) (princ "+"))
             (princ cff)))
    (unless (= pwr 0) (princ "x"))
    (cond ((> pwr 1) (princ (number-superscript pwr))))))

(defun print-polynomial (p)
  (print-monomial (car p))
  (mapcar (lambda (e) (print-monomial e :print-sign t))
          (cdr p)))
