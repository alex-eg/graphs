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

(defun print-matrix (m)
  (let ((row-length (array-dimension m 1)))
    (loop :for i :below (array-total-size m)
       :do (if (zerop (mod i row-length))
               (terpri)
               (princ #\Space))
       (princ (row-major-aref m i))))
    m)

;;;  Matrix (and vector, being 1 x n matrixes) operations

(defun set-edge (m n1 n2 w)
"Sets edge of weighted directed graph with multiedges"
  (push w (aref m n1 n2)))

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

;;; Queue

(defun make-queue (&key (initial-contents))
  initial-contents)

(defun enqueue (q e)
  (append e q))

(defun dequeue (q)
  ())

;;; Search matrix fron position S to all other nodes

(defun breadth-search (s m)
  (let ((nm (make-array (array-dimensions m)
                        :initial-element -1))
        (queue (make-queue :initial-contents (list s))))
    ))


;;; Test

(defun process-case ()
  (let ((m (read-matrix)))
    ))

(let ((num-cases (read)))
  (loop :repeat num-cases
     :do (process-case)))
