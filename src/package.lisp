(cl:in-package :cl)
(defpackage :ceigen-lite (:use))

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defvar ceigen-lite::*src-dir* (asdf:component-pathname (asdf:find-system "ceigen-lite")))
  (defun ceigen-lite::ceigen-lite-c-to-lisp (name)
    (autowrap:default-c-to-lisp
     (cl:if (and (< 11 (length name))
                 (string= "CEIGEN_LITE"
                          (subseq name 0 11)))
            (subseq name 12)
            name)))
  (defun ceigen-lite::starts-with-ceigen-lite (name)
    (alexandria:starts-with-subseq "CEIGEN_LITE_" (symbol-name name)))
  (defun ceigen-lite::does-not-start-with-ceigen-lite (name)
    (not (alexandria:starts-with-subseq "CEIGEN_LITE_" (symbol-name name)))))
