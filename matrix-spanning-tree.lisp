;;; Test for regular and minimal spanning trees count

(defun test (str)
  (let* ((a (cl-math.io:read-weighted-graph-from-string/contest str))
         (mst (cl-math.graph:mst-count-pieper a))
         (st (cl-math.graph:st-count-kirchhoff a)))
    (princ (/ (round mst) (round st)))))
