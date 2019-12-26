;;; Weighted directional graph MST library
;;; Multiple edges are allowed, but loops are not


;;; Supplimentary functions

;;; Matrix I/O
;;; Matrix input is defined by following format:
;;; v e
;;; v1 v2 w1
;;; ...
;;; where v is number of vertices, e - number of edges,
;;; followed by e lines, where v1 and v2 - vertices numbers,
;;; w1 - weight edge

;;; Setters and getters



;;;  Matrix (and vector, being 1 x n matrixes) operations


;;; Various functions for counting MST in a graph
;;; We introduce polynomial matrix, which contais polynomials of arbitrary weight

;;; Kruskal's algorithm related functions

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


;; Polynomial is list of lists, where each member is a pair (power . coefficient)
;; Sorted by power

;;; Test function

(defun test (str)
  (let* ((a (read-matrix-from-string str))
         (mst (mst-count a))
         (st (st-count a)))
    (princ (/ (round mst) (round st)))))
