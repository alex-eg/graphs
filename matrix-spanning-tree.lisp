;;; Weighted directional graph MST library
;;; Multiple edges are allowed, but loops are not


;;; Supplimentary functions

(defun println (obj &optional name)
  "It's like print"
  (if name
      (progn (princ name)
             (princ " ")))
  (princ obj)
  (terpri)
  obj)

(defun string-to-list (string)
  "Read from string as reader would do, and returning it as list.
Example: (string-to-list \"1 2 3 wut quux\") -> '(1 2 3 wut quux)"
    (with-input-from-string (s string)
        (loop
          :for c := (read s nil)
          :while c :collect c)))

(defun make-square-matrix (n)
  (make-array (list n n)
              :initial-element nil))

(defun make-random-square-matrix (n)
  (let ((rm (make-square-matrix n)))
    (loop :for i :below (* n n) :do
       (setf (row-major-aref rm i)
             (let ((r (random 4)))
               (if (= r 0) nil (list r)))))
    rm))

(defun make-random-matrix (n m)
  (let ((rm (make-array (list n m))))
    (loop :for i :below (* n m) :do
       (setf (row-major-aref rm i)
             (random 4)))
    rm))

(defun v-len (v)
  (sqrt (* 1.0d0 ; without double precision results become unstable...
           (loop :for i :below (array-total-size v)
              :sum (let ((a (aref v i)))
                     (* a a))))))

(defun count-population (vec)
  "Counts total length of all vector elements.
Used in counting vertex degree"
  (reduce (lambda (pop e)
            (+ pop (length e)))
          vec
          :initial-value 0))

;;; Matrix I/O
;;; Matrix input is defined by following format:
;;; v e
;;; v1 v2 w1
;;; ...
;;; where v is number of vertices, e - number of edges,
;;; followed by e lines, where v1 and v2 - vertices numbers,
;;; w1 - weight edge

(defun read-matrix-from-string (string)
  (with-input-from-string (*standard-input* string)
    (let* ((nodes (read))
           (edges (read))
           (m (make-square-matrix nodes)))
      (loop :repeat edges
         :do (apply #'set-edge
                    (append (list m)
                            (let ((s (string-to-list (read-line))))
                              (setf (car s) (1- (car s)))
                              (setf (cadr s) (1- (cadr s)))
                              s))))
      m)))

(defun read-matrix ()
  "Reads matrix from *standard-input*"
  (let* ((nodes (read))
         (edges (read))
         (m (make-square-matrix nodes)))
    (loop :repeat edges
       :do (apply #'set-edge
                  (append (list m)
                          (let ((s (string-to-list (read-line))))
                            (setf (car s) (1- (car s)))
                            (setf (cadr s) (1- (cadr s)))
                            s))))
    m))

(defun print-matrix (m &key (print-function #'princ))
  (let ((row-length (array-dimension m 1)))
    (loop :for i :below (array-total-size m)
       :do (if (zerop (mod i row-length))
               (terpri)
               (princ "	"))
       (funcall print-function (row-major-aref m i))))
  m)

;;; Setters and getters

(defun set-edge (m n1 n2 w)
  "Sets edge of weighted directed graph with multiedges"
  (push w (aref m n1 n2))
  (push w (aref m n2 n1)))

(defun set-edge-num (m n1 n2 w)
  "Sets edge of weighted directed graph without multiedges"
  (setf (aref m n1 n2) w)
  (setf (aref m n2 n1) w))

(defun row (m row-num)
  (let ((row-len (array-dimension m 1)))
    (make-array row-len
                :displaced-to m
                :displaced-index-offset (* row-num
                                           row-len))))

(defun column (m col-num)
  (let* ((row-len (array-dimension m 1))
         (col-len (array-dimension m 0))
         (col
          (make-array col-len)))
    (loop :for i :below col-len :do
       (setf (aref col i)
             (row-major-aref m (+ (* row-len i)
                                  col-num))))
    col))

(defun set-colunm (m col-num new-col)
  (let ((n (array-dimension m 0)))
    (loop :for i :below n :do
       (setf (row-major-aref m (+ (* n i)
                                  col-num))
             (aref new-col i)))))

;;;  Matrix (and vector, being 1 x n matrixes) operations

(defmacro defm (name operator)
  "Binary matrix operations definition macro.
Works for weighted directed graphs without multiple edges"
  `(defun ,name (m1 &rest rest)
     (if (null rest) m1
         (let ((m (make-array (array-dimensions m1))))
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

(defun matrix-map (fun m)
  "Maps fun to all m values"
  (let ((nm (make-array (array-dimensions m))))
    (loop :for i :below (array-total-size m)
       :do (setf (row-major-aref nm i)
                 (funcall fun (row-major-aref m i))))
    nm))

(defun m.* (m &rest nums)
  (matrix-map (lambda (e)
                (apply #'* (cons e nums)))
              m))

(defun m./ (m &rest nums)
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

(defun transpose (m1)
  (let ((tm (make-array (nreverse (array-dimensions m1))))
        (n (array-dimension m1 0))
        (m (array-dimension m1 1)))
    (loop :for j :below m
       :do (let ((col (column m1 j)))
             (loop :for i :below n
                :do (setf (aref tm j i) (aref col i)))))
    tm))

(defun number-of-edges (m)
  (loop :for i :below (array-dimension m 0)
     :sum (count-population (row m i))))

(defun vertex-multiplicity (m u v)
  "Number of edges joining u and v, i.e. edge multiplicity"
  (length (aref m u v)))

(defun vertex-degree (m v)
  "Number of edges incident to v"
  (count-population
   (row m v)))

(defun degree-matrix (m)
  (let ((lm (make-array (array-dimensions m)
                        :initial-element 0))
        (n (array-dimension m 1)))
    (loop :for i :below n
       :do (set-edge-num lm i i (vertex-degree m i)))
    lm))

(defun adjacency-matrix (m)
  (let ((lm (make-array (array-dimensions m)
                        :initial-element 0))
        (n (array-total-size m))
        (k (array-dimension m 0)))
    (loop :for i :below n
       :do (setf (row-major-aref lm i)
                 (vertex-multiplicity m (mod i k)
                                      (floor i k))))
    lm))

(defun incidence-matrix (m)
  "Remember that only trimmed of multiedges matrixes are good
for this function"
  (let* ((vertex-num (array-dimension m 0))
         (im (make-array (list vertex-num
                               (number-of-edges m))
                         :initial-element 0))
         (edge-num 0))
    (loop
       :for i :below (array-total-size m)
       :do (let ((current-edge (row-major-aref m i)))
             (if (not (null current-edge))
                 (let ((vertex-from (floor i vertex-num))
                       (vertex-to (mod i vertex-num)))
                   (setf (aref im vertex-from edge-num) -1)
                   (setf (aref im vertex-to edge-num) 1)
                   (incf edge-num)))))
    im))

(defun laplace-matrix (m)
  (m- (degree-matrix m)
      (adjacency-matrix m)))

(defun cofactor (mat row col)
  (let* ((n (array-dimension mat 0))
         (m (array-dimension mat 1))
         (nn (if (> row n)
                 n (1- n)))
         (mm (if (> col m)
                 m (1- m)))
         (mc (make-array (list nn mm))))
    (loop :for i :below nn
       :do (loop
              :for j :below mm
              :do (setf (aref mc i j)
                        (let ((ii (if (>= i row) (1+ i) i))
                              (jj (if (>= j col) (1+ j) j)))
                          (copy-tree (aref mat ii jj))))))
    mc))

(defun dot (u v)
  "u and v dot product"
  (loop :for i :below (array-total-size u)
      :sum (* (aref u i) (aref v i))))

(defun projection (a e)
  "projection of a to e"
  (m.* e (/ (dot a e) (dot e e))))

(defun qr-determinant (m)
  "Find m determinant using qr-decomposition technique"
  (let ((u (make-array (list (array-dimension m 0)))))
    (reduce
     #'*
     (loop :for i :below (array-dimension m 0)
        :collect
        (let* ((ai (column m i))
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

(defun trim-to-min (m)
  "Removes all multi-edges, leaving only minimum value in them"
  (let* ((n (array-dimension m 0))
         (nm (make-array (array-dimensions m)))
         (min-edges 1))
    (loop
       :for i :below n
       :do
       (loop
          :for j :from i :below n
          :do (let*
                  ((edge (aref m i j)))
                (setf (aref nm i j)
                      (setf (aref nm j i)
                            (if (not (null edge))
                                (let ((min (apply #'min edge)))
                                  (setf min-edges (* min-edges
                                                     (count min edge)))
                                  (remove-duplicates
                                   (remove-if-not (lambda (e) (= e min))
                                                  edge)))))))))
    (values nm min-edges)))

;;; Various functions for counting MST in a graph
;;; We introduce polynomial matrix, which contais polynomes of arbitrary weight

;;; Kruskal's algorithm related functions

(defstruct weighted-edge
  "Edge with weight and starting and ending vertices"
  (w 0 :type integer)
  (from 0 :type integer)
  (to 0 :type integer))

(defun weight-edges-to-matrix (we-list)
  "Build an adjacency matrix from weight edges list"
  (let* ((from (mapcar #'weighted-edge-from we-list))
         (to (mapcar #'weighted-edge-to we-list))
         (vert-num (1+ (apply #'max (nconc from to))))
         (m (make-array (list vert-num vert-num))))
    (dolist (e we-list)
      (setf (aref m (weighted-edge-from e)
                  (weighted-edge-to e))
            (weighted-edge-w e)))
    m))

(defun edges-vertices-list (m)
  "Having graph matrix, convert it to list of weighted edges"
  (let (edge-vert-list)
    (loop :for i :below (array-total-size m)
       :do (let* ((weight (row-major-aref m i))
                  (n (array-dimension m 0))
                  (from (floor i n))
                  (to (mod i n)))
             (when weight
               (push (make-weighted-edge :w (car weight)
                                         :from from
                                         :to to)
                     edge-vert-list))))
    edge-vert-list))

(defun subset-find (set v)
  "Finds a subset containing and element v"
  (find v set :test (lambda (v set)
                      (find v set))))

(defun subset-remove (set subset)
  (delete subset set :test #'eq))

(defun subset-merge (set subset-1 subset-2)
  (nconc subset-1 subset-2)
  (subset-remove set subset-2))

(defun find-mst (m)
  "Find mst using Kruskal's algorithm"
  (let* ((edges (sort (edges-vertices-list (trim-to-min m))
                      (lambda (a b) (< (weighted-edge-w a)
                                       (weighted-edge-w b)))))
         (verts-num (array-dimension m 0))
         (subsets (loop :for i :below verts-num :collect (list i)))
         mst)
    (dolist (e edges)
      (let* ((from (weighted-edge-from e))
             (to (weighted-edge-to e))
             (from-set (subset-find subsets from))
             (to-set (subset-find subsets to)))
        (unless (eq from-set to-set)    ; Otherwise they would form a cycle
          (setf subsets (subset-merge subsets from-set to-set))
          (push e mst))))
    mst))


(defun find-all (item list &key (key #'identity) (test #'=))
  (let (res)
    (dolist (e list res)
      (if (funcall test (funcall key e) item)
          (push e res)))))

(defun remove-all (item-list list)
  (reduce (lambda (a e) (remove e a))
          item-list :initial-value list))

(defun find-item (tree num)
  "Tree is classic: car is value, cdr is list with all children"
  (if (= num (car tree))
      tree
      (car
       (remove-if #'null
                  (mapcar (lambda (subtree) (find-item subtree num))
                          (cdr tree))))))

(defun build-tree-from-edges (weighted-edges)
  (let* ((from (mapcar #'weighted-edge-from weighted-edges))
         (to (mapcar #'weighted-edge-to weighted-edges))
         (root (apply #'max (append from to)))
         (edges (copy-list weighted-edges)))
    (labels ((find-and-delete-all-children (node)
               (let* ((tos (find-all node edges :key #'weighted-edge-to))
                      (froms (find-all node edges :key #'weighted-edge-from))
                      (children (append (mapcar #'weighted-edge-from tos)
                                        (mapcar #'weighted-edge-to froms))))
                 (setf edges (remove-all (append tos froms) edges))
                 children))
             (build-tree (node)
               (let ((children (find-and-delete-all-children  (car node))))
                 (setf (cdr node) (mapcar #'list children))
                 (mapcar #'build-tree (cdr node)))))
      (let ((tree (list root)))
        (build-tree tree)
        tree))))

(defun find-leaf (tree)
  (let ((children (cdr tree)))
    (if (null children)
        tree
        (if (null (cdr (car children)))
            tree
            (find-leaf (car children))))))

;;; Pieper algorithm for counting MST in graph

(defun a~-matrix (m n)
  "Matrix m without the n-th row"
  (let* ((col-num (array-dimension m 1)))
    (cofactor m n (1+ col-num))))

;; Polynome is list of lists, where each member is a pair (power . coefficient)
;; Sorted by power

(defun number-superscript (number)
  (coerce (nreverse (%number-superscript number)) 'string))

(defun %number-superscript (number)
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
              (%number-superscript (floor number 10))))))

(defun print-monome (m &key print-sign)
  (let ((pwr (pwr m))
        (cff (cff m)))
    (cond ((= cff 0) (princ 0) (return-from print-monome))
          ((and (= cff 1) print-sign) (princ "+"))
          ((= cff -1) (princ "-"))
          (t (if (and (> cff 1) print-sign) (princ "+"))
           (princ cff)))
    (unless (= pwr 0) (princ "x"))
    (cond ((> pwr 1) (princ (number-superscript pwr))))))

(defun print-polynome (p)
  (print-monome (car p))
  (mapcar (lambda (e) (print-monome e :print-sign t))
          (cdr p)))

(defun make-monome (power coefficient)
  (cons power coefficient))

(defun pwr (a)
"Monome power"
  (car a))

(defun cff (a)
"Monome coefficient"
  (cdr a))

(defun polynome+ (a b)
  (let* ((sum (%polynome+ a b)))
    (setf sum (remove-if #'zerop sum :key #'cff))
    (if (null sum)
        (list (make-monome 0 0))
        sum)))

(defun %polynome+ (a b)
  (cond ((null b) a)
        ((null a) b)
        ((> (pwr (car a))
            (pwr (car b)))
         (cons (car a) (%polynome+ (cdr a) b)))
        ((< (pwr (car a))
            (pwr (car b)))
         (cons (car b) (%polynome+ a (cdr b))))
        ((= (pwr (car a))
            (pwr (car b)))
         (cons (make-monome (pwr (car a))
                            (+ (cff (car a))
                               (cff (car b))))
               (%polynome+ (cdr a) (cdr b))))))

(defun polynome/ (polynome monome)
  (let ((polynome (copy-tree polynome))
        (monome (copy-tree monome)))
    (mapcar (lambda (mm)
              (if (= (cff mm) 0)
                  mm
                  (progn
                    (rplaca mm (- (pwr mm) (pwr monome)))
                    (rplacd mm (/ (cff mm) (cff monome))))))
            polynome)))

(defun polynome-min-power (polynome)
  (apply #'min (mapcar #'car (mapcar #'last polynome))))

(defun monome* (a num)
  (make-monome (pwr a)
                 (* num (cff a))))

(defun vertex-degree-polynome (m u)
  "Sum of all weigth polynomes of the vertex"
  (let* ((edges-starting (row m u))
         (polynomes (mapcar (lambda (a) (make-monome a 1))
                            (reduce #'append
                                 (remove-if #'null edges-starting)))))
    (if polynomes
        (sort (reduce (lambda (acc p)
                        (polynome+ acc (list p)))
                       (cdr polynomes)
                       :initial-value (list (car polynomes)))
              #'> :key #'car)
        (list (make-monome 0 0)))))

(defun edge-list-weight-polynome (m u v)
  (let ((edge-list (mapcar (lambda (a) (make-monome a -1))
                           (aref m u v))))
    (if edge-list
        (sort (reduce (lambda (acc p)
                        (polynome+ acc (list p)))
                      (cdr edge-list)
                      :initial-value (list (car edge-list)))
              #'> :key #'car)
        (list (make-monome 0 0)))))

(defun weighted-adjacency-matrix (m)
  (let ((lm (make-array (array-dimensions m)))
        (n (array-dimension m 0)))
    (loop :for i :below n
       :do (loop :for j :below n
              :do
              (setf (aref lm i j)
                    (if (= i j)
                        (vertex-degree-polynome m i)
                        (edge-list-weight-polynome m i j)))))
    lm))

(defun trace-weighted-mst (wm mst)
  "The key function to Pieper algorithm. We traverse MST and modify weighted
matrix by linear combination of columns"
  (let* ((n (array-dimension wm 0))
         (co-wm (cofactor wm n n)))
    (do ((f (find-leaf mst) (find-leaf mst)))
        ((= 1 (length mst)))
      (let* ((parent (car f))
             (leaf (caadr f)))
        (unless (or (= parent (1- n))
                    (= leaf (1- n)))
          (let ((col-1 (column co-wm parent))
                (col-2 (column co-wm leaf)))
;            (format t "~%~%~A->~A~%" leaf parent)
            (set-colunm co-wm parent
                        (map 'vector (lambda (e)
                                       (if (null e)
                                           (list (make-monome 0 0))
                                           e))
                             (map 'vector #'polynome+ col-1 col-2)))
 ;           (print-matrix co-wm :print-function #'print-polynome)
            ))
        (setf (cdr f) (cddr f))))
    co-wm))

(defun factor-out (m)
  "Factors out monomials and returns matrix with value x = 0"
  (let ((n (array-dimension m 0))
        (nm (make-array (array-dimensions m))))
    (loop :for j :below n :do
       (let* ((col (column m j))
              (min-pow (apply #'min
                              (remove 0
                                      (map 'list #'polynome-min-power col)))))
         (loop :for i :below n :do
            (setf (aref nm i j)
                  (let ((min (car (last
                                   (polynome/ (aref col i)
                                              (cons min-pow 1))))))
                    (if (= 0 (pwr min))
                        (cff min)
                        0))))))
    nm))

(defun mst-count (m)
  (qr-determinant
   (factor-out
    (trace-weighted-mst
     (weighted-adjacency-matrix m)
     (build-tree-from-edges (find-mst m))))))

(defun st-count (m)
  (qr-determinant (cofactor (laplace-matrix m) 0 0)))


;;; Test function

(defun test ()
  (princ (/ 1 (span-tree-number (read-matrix-from-string "4 4
1 2 1
2 3 4
3 4 3
4 1 2")))))
