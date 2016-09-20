;;; Weighted directional graph MST library
;;; Multiple edges are allowed, but loops are not


;;; Supplimentary functions

(defun string-to-list (string)
  "Read from string as reader would do, and returning it as list.
Example: (string-to-list \"1 2 3 wut quux\") -> '(1 2 3 wut quux)"
    (with-input-from-string (s string)
        (loop
          :for c := (read s nil)
          :while c :collect c)))

;; Matrix I/O
;;; Matrix input is defined by following format:
;;; v e
;;; v1 v2 w1
;;; ...
;;; where v is number of vertices, e - number of edges,
;;; followed by e lines, where v1 and v2 - vertices numbers,
;;; w1 - weight edge

(defun make-square-matrix (n)
  (make-array (list n n) :initial-element nil))

(defun read-matrix ()
  "Reads matrix from *standard-input*"
  (let* ((nodes (read))
         (edges (read))
         (m (make-square-matrix nodes)))
    (loop :repeat edges
       :do (let ((s (string-to-list (read-line))))
             (setf (aref m
                         (1- (car s))
                         (1- (cadr s)))
                   6
                   (aref m
                         (1- (cadr s))
                         (1- (car s)))
                   6)))
    m))

;;;  Matrix (and vector, being 1 x n matrixes) operations

(defun row (m row-num)
  (let ((row-len (array-dimension m 1)))
    (make-array row-len
                :displaced-to m
                :displaced-index-offset (* row-num
                                           row-len))))

;;; Queue

(defun make-queue ()
  (cons nil nil))

(defun enqueue (obj q)
  (if (null (car q))
      (setf (cdr q) (setf (car q) (list obj)))
      (setf (cdr (cdr q)) (list obj)
            (cdr q) (cdr (cdr q))))
  (car q))

(defun queue-nullp (q)
  (null (car q)))

(defun dequeue (q)
  (pop (car q)))

;;; Search matrix fron position S to all other nodes

(defun breadth-search (m s)
  (let ((nm (make-array (array-dimension m 0) :initial-element -1))
        (q (make-queue)))
    (enqueue (cons s 0) q)
    (loop :while (not (queue-nullp q))
       :do (let* ((node-dist (dequeue q))
                  (current-node (car node-dist))
                  (current-distance (cdr node-dist))
                  (adjacent (row m current-node)))
             (loop :for i :below (array-dimension m 0)
                :do (let ((edge (aref adjacent i)))
                      (if edge
                          (if (and (= -1 (aref nm i))
                                   (not (= i s)))
                              (progn
                                (setf (aref nm i) (+ 6 current-distance))
                                (enqueue (cons i (+ 6 current-distance)) q))))))))
    (loop :for i :below (array-dimension nm 0)
       :when (not (= i s)) :collect (aref nm i))))


;;; Test

(defun process-case ()
  (let ((m (read-matrix))
        (start (1- (read))))
    (format t "~{~A~^ ~}~%" (breadth-search m start))))

(let ((num-cases (read)))
  (loop :repeat num-cases
     :do (process-case)))
