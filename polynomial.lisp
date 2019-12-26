(defpackage :cl-math.polynomial
  (:use :cl)
  (:export
   :make-polynomial
   ))

(in-package :cl-math.polynomial)

;;; Polynomial is list of lists, where each member is a monomial,
;;; which is a cons of following structure: (power . coefficient)
;;; Monomials in a polynomial are sorted by power

(defun make-monomial (power coefficient)
  "Make "
  (cons power coefficient))

(defun pwr (a)
  "Monomial power"
  (car a))

(defun cff (a)
  "Monomial coefficient"
  (cdr a))

(defun make-polynomial (&rest monomials)
  (apply #'list (sort monomials #'> :key #'pwr)))

(defun polynomial-min-power (polynomial)
  (apply #'min (mapcar #'car (mapcar #'last polynomial))))

(defun monomial* (a num)
  "Multiply monomial by number"
  (make-monomial (pwr a)
                 (* num (cff a))))

(defun polynomial+ (a b)
  (cond
    ((and (null a)
          (null b))
     (list (make-monomial 0 0)))
    ((null b) a)
    ((null a) b)
    ((> (pwr (car a))
        (pwr (car b)))
     (cons (car a) (polynomial+ (cdr a) b)))
    ((< (pwr (car a))
        (pwr (car b)))
     (cons (car b) (polynomial+ a (cdr b))))
    ((= (pwr (car a))
        (pwr (car b)))
     (let ((new-monomial
             (make-monomial (pwr (car a))
                            (+ (cff (car a))
                               (cff (car b))))))
       (if (= 0 (cff new-monomial))
         (polynomial+ (cdr a) (cdr b))
         (cons new-monomial
               (polynomial+ (cdr a) (cdr b))))))))

(defun polynomial/ (polynomial monomial)
  (let ((polynomial (copy-tree polynomial))
        (monomial (copy-tree monomial)))
    (mapcar (lambda (mm)
              (if (= (cff mm) 0)
                mm
                (progn
                  (rplaca mm (- (pwr mm) (pwr monomial)))
                  (rplacd mm (/ (cff mm) (cff monomial))))))
            polynomial)))
