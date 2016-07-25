;; Enter your code here. Read input from STDIN. Print output to STDOUT

(defun string-to-list (string)
    (with-input-from-string (s string)
        (loop
          :for c := (read s nil)
          :while c :collect c)))

(defun make-square-matrix (n)
  (make-array (list n n) :element-type 'integer))

(defun print-matrix (m)
  (let ((row-length (array-dimension m 1)))
    (loop :for i :below (array-total-size m)
       :do (if (zerop (mod i row-length))
               (terpri)
               (princ #\Space))
       (princ (row-major-aref m i))))
    m)

(defun set-edge (m n1 n2 w)
  (setf (aref m n1 n2) w)
  (setf (aref m n2 n1) w))

(defun get-row (m row-num)
  (let ((row-len (array-dimension m 1)))
    (make-array row-len
                :element-type 'integer
                :displaced-to m
                :displaced-index-offset (* row-num
                                           row-len))))

(defun get-column (m col-num)
  (let* ((row-len (array-dimension m 1))
         (col-len (array-dimension m 0))
         (col
          (make-array row-len
                      :element-type 'integer)))
    (loop :for i :below col-len :do
       (setf (aref col i)
             (row-major-aref m (+ (* row-len i)
                                  col-num))))
    col))

(defun v-len (v)
  (sqrt (loop :for i :below (array-total-size v)
           :sum (expt (aref v i) 2))))

(defun count-population (vec)
  (loop :for v :across vec :counting (> v 0)))

(defmacro defm (name operator)
  `(defun ,name (m1 &rest rest)
     (if (null rest) m1

         (let ((m (make-array (array-dimensions m1)
                              :element-type 'integer)))
           (loop :for i :below (array-total-size m)
              :do (setf (row-major-aref m i)
                        (apply ,operator
                               (cons
                                (row-major-aref m1 i)
                                (mapcar (lambda (m) (row-major-aref m i))
                                        rest)))))
           m))))

(defm m+ #'+)
(defm m- #'-)

(defmacro defm. (name op)
  `(defun ,name (m &rest nums)
     (if (null nums) m
         (let ((nm (make-array (array-dimensions m)
                               :element-type 'integer)))
           (loop :for i :below (array-total-size m)
              :do (setf (row-major-aref nm i)
                        (apply ,op
                               (cons
                                (row-major-aref m i)
                                nums))))
           nm))))

(defm. m.* #'*)
(defm. m./ #'/)

(defun degree-matrix (m)
  (let ((lm (make-array (array-dimensions m)
                        :element-type 'integer
                        :initial-element 0))
        (n (array-dimension m 1)))
    (loop :for i :below n
       :do (set-edge lm i i (count-population
                             (get-row m i))))
    lm))

(defun adjacency-matrix (m)
  (let ((am (make-array (array-dimensions m)
                        :element-type 'integer
                        :initial-element 0))
        (n (array-total-size m)))
    (loop :for i :below n
       :do (setf (row-major-aref am i)
                 (if (zerop (row-major-aref m i))
                     0 1)))
    am))

(defun laplace-matrix (m)
  (m- (degree-matrix m)
      (adjacency-matrix m)))

(defun cofactor (mat row col)
  (let* ((mc (make-array (mapcar #'1- (array-dimensions mat))
                         :element-type 'integer))
         (n (array-dimension mc 0))
         (m (array-dimension mc 1)))
    (loop :for i :below n
       :do (loop
              :for j :below m
              :do (setf (aref mc i j)
                        (let ((ii (if (>= i row) (1+ i) i))
                              (jj (if (>= j col) (1+ j) j)))
                          (aref mat ii jj)))))
    mc))

(defun determinant (m)
  (cond ((equal (array-dimensions m) '(2 2))
         (- (* (aref m 0 0) (aref m 1 1))
            (* (aref m 0 1) (aref m 1 0))))
        ((equal (array-dimensions m) '(1 1))
         (aref m 0 0))
        (t (loop :for i :below (array-dimension m 0)
              :sum (* (if (evenp i) 1 -1)
                      (aref m 0 i)
                      (determinant (cofactor m 0 i)))))))

(defun dot (u v)
  (loop :for i :below (array-total-size u)
      :sum (* (aref u i) (aref v i))))

(defun projection (a e)
  (m.* e (/ (dot a e) (dot e e))))

(defun println (obj &optional name)
  (if name
      (progn (princ name)
             (princ " ")))
  (princ obj)
  (terpri)
  obj)

(defun qr-determinant (m)
  (let ((u (make-array (list (array-dimension m 0)))))
    (reduce #'*
            (loop :for i :below (array-dimension m 0)
               :collect (let* ((ai (get-column m i))
                               (ui
                                (setf (aref u i)
                                      (apply
                                       #'m-
                                       (cons
                                        ai
                                        (loop :for j :from 1 :to i
                                           :collect
                                           (projection ai (aref u (1- j))))))))
                               (ei (m./ ui (v-len ui))))
                          (dot ei ai)))
            :initial-value 1)))

(defun min-edges-num (m)
  (loop :for i :below (array-dimension m 0)
     :collect ()))

(defun span-tree-number (m)
  (determinant (cofactor (laplace-matrix m) 3 3)))

(defun random-matrix (n)
  (let ((rm (make-square-matrix n)))
    (loop :for i :below (* n n) :do
       (setf (row-major-aref rm i) (random 4)))
    rm))

(defun read-matrix-from-string (string)
  (with-input-from-string (*standard-input* string)
    (let* ((nodes (read))
           (edges (read))
           (m (make-square-matrix nodes)))
      (loop :repeat edges
         :do (apply #'set-edge (append (list m)
                                       (let ((s (string-to-list (read-line))))
                                         (setf (car s) (1- (car s)))
                                         (setf (cadr s) (1- (cadr s)))
                                         s))))
      m)))

(defun test ()
  (princ (/ 1 (span-tree-number (read-matrix-from-string "4 4
1 2 1
2 3 4
3 4 3
4 1 2")))))
