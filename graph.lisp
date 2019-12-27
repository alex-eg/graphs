(defpackage :cl-math.graph
  (:use
   :cl
   :cl-math.polynomial
   :cl-math.matrix
   :cl-math.list
   :cl-math.set)

  (:export
   :set-edge
   :set-edge-num
   :number-of-edges
   :vertex-multiplicity
   :vertex-degree

   ;; Graph matrices
   :degree-matrix
   :adjacency-matrix
   :incidence-matrix
   :laplace-matrix

   ;; Spanning tree algorithms
   :mst-find-kruskal
   :st-count-kirchhoff
   :mst-count-pieper))

(in-package :cl-math.graph)

(defun set-edge (m n1 n2 w)
  "Set edge of weighted directed graph with multiedges"
  (push w (aref m n1 n2))
  (push w (aref m n2 n1)))

(defun set-edge-num (m n1 n2 w)
  "Set edge of weighted directed graph without multiedges"
  (setf (aref m n1 n2) w)
  (setf (aref m n2 n1) w))

(defun number-of-edges (m)
  "Get number of edges in the graph"
  (loop :for i :below (array-dimension m 0)
     :sum (count-population (row m i))))

(defun vertex-multiplicity (m u v)
  "Get number of edges joining u and v, i.e. edge multiplicity"
  (length (aref m u v)))

(defun vertex-degree (m v)
  "Number of edges incident to v"
  (count-population
   (row m v)))

(defun degree-matrix (m)
  (let ((lm (make-array (array-dimensions m)
                        :initial-element 0))
        (n (array-dimension m 1)))
    (loop
      :for i :below n
      :do (set-edge-num lm i i (vertex-degree m i)))
    lm))

(defun adjacency-matrix (m)
  (let ((lm (make-array (array-dimensions m)
                        :initial-element 0))
        (n (array-total-size m))
        (k (array-dimension m 0)))
    (loop
      :for i :below n
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

(defun laplace-matrix (g)
  "Compute laplacian matrix of the graph"
  (m- (degree-matrix g)
      (adjacency-matrix g)))

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
           :do
              (let* ((edge (aref m i j)))
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

;;; Another form of graph representation: using list of weighted edges
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

(defun mst-find-kruskal (m)
  "Find (a certain) MST using Kruskal's algorithm"
  (let* ((edges (sort (edges-vertices-list (trim-to-min m))
                      (lambda (a b) (< (weighted-edge-w a)
                                       (weighted-edge-w b)))))
         (verts-num (array-dimension m 0))
         (subsets (loop :for i :below verts-num :collect (list i)))
         (mst (list)))
    (dolist (e edges)
      (let* ((from (weighted-edge-from e))
             (to (weighted-edge-to e))
             (from-set (subset-find subsets from))
             (to-set (subset-find subsets to)))
        (unless (eq from-set to-set) ; Otherwise they would form a cycle
          (setf subsets (subset-merge subsets from-set to-set))
          (push e mst))))
    mst))

;;; Pieper algorithm for counting MST in graph
(defun vertex-degree-polynomial (m u)
  "Sum of all weigth polynomials of the vertex"
  (let* ((edges-starting (row m u))
         (polynomials (mapcar (lambda (a) (make-monomial a 1))
                              (reduce #'append
                                      (remove-if #'null edges-starting)))))
    (if polynomials
      (sort (reduce (lambda (acc p)
                      (polynomial+ acc (list p)))
                    (cdr polynomials)
                    :initial-value (list (car polynomials)))
            #'> :key #'car)
      (list (make-monomial 0 0)))))

(defun edge-list-weight-polynomial (m u v)
  (let ((edge-list (mapcar (lambda (a) (make-monomial a -1))
                           (aref m u v))))
    (if edge-list
        (sort (reduce (lambda (acc p)
                        (polynomial+ acc (list p)))
                      (cdr edge-list)
                      :initial-value (list (car edge-list)))
              #'> :key #'car)
        (list (make-monomial 0 0)))))

(defun weighted-adjacency-matrix (m)
  (let ((lm (make-array (array-dimensions m)))
        (n (array-dimension m 0)))
    (loop
      :for i :below n
      :do (loop :for j :below n
                :do
                   (setf (aref lm i j)
                         (if (= i j)
                           (vertex-degree-polynomial m i)
                           (edge-list-weight-polynomial m i j)))))
    lm))

(defun trace-weighted-mst (wm mst)
  "The key function to Pieper algorithm. We traverse MST and modify weighted
matrix by linear combination of columns"
  (let* ((n (array-dimension wm 0))
         (co-wm (cofactor-matrix wm n n)))
    (do ((f (find-leaf mst) (find-leaf mst)))
        ((= 1 (length mst)))
      (let* ((parent (car f))
             (leaf (caadr f)))
        (unless (or (= parent (1- n))
                    (= leaf (1- n)))
          (let ((col-1 (column co-wm parent))
                (col-2 (column co-wm leaf)))
            (set-column co-wm parent
                        (map 'vector
                             (lambda (e)
                               (if (null e)
                                 (list (make-monomial 0 0))
                                 e))
                             (map 'vector #'polynomial+ col-1 col-2)))))
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
                                     (map 'list #'polynomial-min-power col)))))
        (loop :for i :below n :do
          (setf (aref nm i j)
                (let ((min (car (last
                                 (polynomial/ (aref col i)
                                              (cons min-pow 1))))))
                  (if (= 0 (pwr min))
                    (cff min)
                    0))))))
    nm))

(defun find-item (tree num)
  "Tree is classic: car is value, cdr is list with all children"
  (if (= num (car tree))
    tree
    (car
     (remove-if #'null
                (mapcar (lambda (subtree) (find-item subtree num))
                        (cdr tree))))))

(defun find-leaf (tree)
  "Find arbitrary leaf in a tree"
  (let ((children (cdr tree)))
    (if (null children)
      tree
      (if (null (cdr (car children)))
        tree
        (find-leaf (car children))))))

(defun build-tree-from-edges (weighted-edges)
  "Supplimentary function to build tree fromm list of graph weighted edges"
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
               (let ((children (find-and-delete-all-children (car node))))
                 (setf (cdr node) (mapcar #'list children))
                 (mapcar #'build-tree (cdr node)))))
      (let ((tree (list root)))
        (build-tree tree)
        tree))))

(defun mst-count-pieper (m)
  "Find number of minimal spanning trees in a weighted graph using Pieper algorithm"
  (gj-determinant
   (factor-out
    (trace-weighted-mst
     (weighted-adjacency-matrix m)
     (build-tree-from-edges
      (mst-find-kruskal m))))))

(defun st-count-kirchhoff (m)
  "Find number of all spanning trees in a graph using Kirchhoff's theorem"
  (gj-determinant (cofactor-matrix (laplace-matrix m) 0 0)))
