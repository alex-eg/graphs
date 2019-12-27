(defpackage :cl-math.polynomial
  (:use :cl)
  (:export
   :make-monomial
   :pwr
   :cff
   :make-polynomial
   :polynomial+
   :polynomial/
   :polynomial-min-power))

(in-package :cl-math.polynomial)

;;; Polynomial is list of lists, where each member is a monomial,
;;; which is a cons of following structure: (power . coefficient)
;;; Monomials in a polynomial are sorted by power

(defun make-monomial (power coefficient)
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

(defun polynomial-zerop (a)
  (and (= 1 (length a))
       (= (cff (car a))
          (pwr (car a))
          0)))

(defun polynomial+ (a b)
  (labels
      ((%+ (a b)
         (cond
           ((polynomial-zerop a) b)
           ((polynomial-zerop b) a)
           ((null b) a)
           ((null a) b)
           ((> (pwr (car a))
               (pwr (car b)))
            (cons (car a) (%+ (cdr a) b)))
           ((< (pwr (car a))
               (pwr (car b)))
            (cons (car b) (%+ a (cdr b))))
           ((= (pwr (car a))
               (pwr (car b)))
            (let ((new-monomial
                    (make-monomial (pwr (car a))
                                   (+ (cff (car a))
                                      (cff (car b))))))
              (if (= 0 (cff new-monomial))
                (%+ (cdr a) (cdr b))
                (cons new-monomial
                      (%+ (cdr a) (cdr b)))))))))
    (let ((sum (%+ a b)))
      (if (null sum)
        (make-polynomial (make-monomial 0 0))
        sum))))

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
