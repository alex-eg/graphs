;;; Also contains vector operations

(defpackage :cl-math.matrix
  (:use :cl)
  (:export
   :make-matrix
   :make-square-matrix
   :make-random-matrix
   :make-random-square-matrix

   :row
   :set-row
   :column
   :set-column

   :m+
   :m-
   :m.*
   :m./
   :m*

   ;; Vector
   :v-len
   :dot
   :projection
   :count-population

   :transpose
   :cofactor-matrix
   :a~-matrix
   :gauss-jordan-eliminate
   :determinant
   :qr-determinant
   :qr-determinant-gram-schmidt
   :gj-determinant))

(in-package :cl-math.matrix)

(defun make-matrix (n m)
  (make-array (list n m)
              :initial-element nil))

(defun make-square-matrix (n)
  (make-matrix n n))

(defun make-random-matrix (n m &key (random-cap 4))
  (let ((rm (make-array (list n m))))
    (loop
      :for i :below (* n m)
      :do (setf (row-major-aref rm i)
                (random random-cap)))
    rm))

(defun make-random-square-matrix (n)
  (make-random-matrix n n))

(defun row (m row-num)
  (let ((row-len (array-dimension m 1)))
    (make-array row-len
                :displaced-to m
                :displaced-index-offset (* row-num
                                           row-len))))
(defun set-row (m row-num new-row)
  "Set matrix row to new value.
New-row length must be greater or equal than matrix 1-st dimension"
  (let ((row-len (array-dimension m 1)))
    (loop
      :for i :from (* row-num row-len)
        :below (* (1+ row-num) row-len)
      :for k :from 0 :do
        (setf (row-major-aref m i) (aref new-row k)))))

(defun column (m col-num)
  (let* ((row-len (array-dimension m 1))
         (col-len (array-dimension m 0))
         (col (make-array col-len)))
    (loop
      :for i :below col-len :do
        (setf (aref col i)
              (row-major-aref m (+ (* row-len i)
                                   col-num))))
    col))

(defun set-column (m col-num new-col)
  "Set matrix columnt to new value.
New-col length must be greater or equal than matrix 2-nd dimension"
  (let ((n (array-dimension m 0)))
    (loop
      :for i :below n :do
        (setf (row-major-aref m (+ (* n i)
                                   col-num))
              (aref new-col i)))))

(defmacro defm (name operator)
  "Binary matrix operations definition macro.
Works for weighted directed graphs without multiple edges"
  `(defun ,name (m1 &rest rest)
     (if (null rest) m1
       (let ((m (make-array (array-dimensions m1))))
         (loop
           :for i :below (array-total-size m)
           :do (setf (row-major-aref m i)
                     (apply ,operator
                            (cons
                             (row-major-aref m1 i)
                             (mapcar (lambda (m) (row-major-aref m i))
                                     rest)))))
         m))))

(defm m+ #'+)
(defm m- #'-)

(defun matrix-map (fun m)
  "Map fun to all elements of matrix m"
  (let ((nm (make-array (array-dimensions m))))
    (loop :for i :below (array-total-size m)
          :do (setf (row-major-aref nm i)
                    (funcall fun (row-major-aref m i))))
    nm))

(defun m.* (m &rest nums)
  "Multiply each element of m by all of nums"
  (matrix-map (lambda (e)
                (apply #'* (cons e nums)))
              m))

(defun m./ (m &rest nums)
  "Divide each element of m by all of nums"
  (matrix-map (lambda (e)
                (apply #'/ (cons e nums)))
              m))

(defun m* (m1 m2)
  "Regular matrix multiplication"
  (assert (equal (array-dimensions m1)
                 (nreverse (array-dimensions m2))))
  (let* ((n (array-dimension m1 0))
         (m (array-dimension m2 1))
         (nm (make-array (list n m))))
    (dotimes (i n)
      (dotimes (j m)
        (setf (aref nm i j)
              (reduce (lambda (acc elems)
                        (+ acc (* (car elems)
                                  (cdr elems))))
                      (map 'list #'cons
                           (row m1 i)
                           (column m2 j))
                      :initial-value 0))))
    nm))

(defun v-len (v)
  "Compute Euclidian norm of a vector"
  (sqrt (* 1.0d0 ; without double precision results become unstable...
           (loop
             :for i :below (array-total-size v)
             :sum (let ((a (aref v i)))
                    (* a a))))))

(defun dot (u v)
  "Compute u and v dot product"
  (loop
    :for i :below (array-total-size u)
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

(defun transpose (mat)
  "Regular matrix transposition"
  (let ((tm (make-array (nreverse (array-dimensions mat))))
        (n (array-dimension mat 0))
        (m (array-dimension mat 1)))
    (loop
      :for j :below m
      :do (let ((col (column mat j)))
            (loop
              :for i :below n
              :do (setf (aref tm j i) (aref col i)))))
    tm))

(defun cofactor-matrix (mat row col)
  "Find cofactor matrix"
  (let* ((n (array-dimension mat 0))
         (m (array-dimension mat 1))
         (nn (if (> row n)
               n (1- n)))
         (mm (if (> col m)
               m (1- m)))
         (mc (make-array (list nn mm))))
    (loop
      :for i :below nn
      :do (loop
            :for j :below mm
            :do (setf (aref mc i j)
                      (let ((ii (if (>= i row) (1+ i) i))
                            (jj (if (>= j col) (1+ j) j)))
                        (copy-tree (aref mat ii jj))))))
    mc))

(defun a~-matrix (m n)
  "Matrix m without the n-th row"
  (let* ((col-num (array-dimension m 1)))
    (cofactor m n (1+ col-num))))

(defun determinant (m)
  "Find determinant using (modified) naive implementation"
  (let ((n (array-dimension m 0)))
    (cond ((= n 1) (row-major-aref m 0))
          ((= n 2)
           (- (* (row-major-aref m 0)
                 (row-major-aref m 3))
              (* (row-major-aref m 1)
                 (row-major-aref m 2))))
          ((> n 2)
           (apply #'+
                  (loop
                    :for i :below n
                    :collect (if (= (row-major-aref m i) 0)
                               0
                               (* (if (oddp i) -1 1)
                                  (row-major-aref m i)
                                  (determinant (cofactor m 0 i))))))))))


(defun compute-u-k (k v u-vector)
  "Gram-Schmidt helper for computing k-th u vector"
  (if (= k 0) v
    (loop
      :for uk := (m- v (projection v (aref u-vector j)))
        :then (m- uk (projection uk (aref u-vector j)))
      :for j :below k
      :finally (return uk))))

(defun qr-determinant-gram-schmidt (m)
  "Find m determinant using QR decomposition technique,
with Modified Gram-Schmidt process"
  (let ((u (make-array (list (array-dimension m 0)))))
    (reduce
     #'*
     (loop
       :for i :below (array-dimension m 0)
       :collect
       (let* ((ai (column m i))
              (ui
                (setf (aref u i)
                      (compute-u-k i ai u)))
              (ei (m./ ui (v-len ui))))
         (dot ei ai)))
     :initial-value 1)))

(defun qr-determinant (m)
  "Find m determinant using QR decomposition technique"
  (let ((u (make-array (list (array-dimension m 0)))))
    (reduce
     #'*
     (loop
       :for i :below (array-dimension m 0)
       :collect
       (let* ((ai (column m i))
              (ui
                (setf (aref u i)
                      (apply
                       #'m-
                       (cons
                        ai
                        (loop
                          :for j :below i
                          :collect
                          (projection ai (aref u j)))))))
              (ei (m./ ui (v-len ui))))
         (dot ei ai)))
     :initial-value 1)))

(defun gauss-jordan-eliminate (m)
  "Coerce matrix to echelon form using Gauss-Jordan elimination"
  (let ((n (array-dimension m 0)))
    (loop
      :for i :below n
      :do
         (let ((row (row m i)))
           (loop
             :for ii :from (1+ i) :below n
             :do
                (let* ((current-row (row m ii))
                       (multiplier (/ (aref current-row i)
                                      (aref row i)))
                       (sub-row (m.* row multiplier))
                       (new-row (m- current-row sub-row)))
                  (set-row m ii new-row)))))
    m))

(defun gj-determinant (m)
  "Find determinant using Gauss-Jordan elimination"
  (let ((gj-m (gauss-jordan-eliminate m))
        (det 1))
    (loop
      :for i :below (array-dimension m 0)
      :do (setf det (* det (aref gj-m i i))))
    det))
