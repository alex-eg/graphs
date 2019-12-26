(defpackage :cl-math.vector
  (:use :cl)
  (:export
   :dot
   :projection
   ))

(defun v-len (v)
  "Compute Euclidian norm of a vector"
  (sqrt (* 1.0d0 ; without double precision results become unstable...
           (loop
             :for i :below (array-total-size v)
             :sum (let ((a (aref v i)))
                    (* a a))))))

(defun dot (u v)
  "Compute u and v dot product"
  (loop :for i :below (array-total-size u)
      :sum (* (aref u i) (aref v i))))

(defun projection (u v)
  "Find projection of u to v"
  (m.* e (/ (dot a e) (dot e e))))

(defun count-population (vec)
  "Counts total length of all vector elements.
Used in counting vertex degree of weighted graph"
  (reduce (lambda (pop e)
            (+ pop (length e)))
          vec
          :initial-value 0))
