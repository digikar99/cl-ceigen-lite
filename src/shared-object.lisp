(cl:in-package :cl)

(eval-when (:compile-toplevel :load-toplevel :execute)
  (do-symbols (s (find-package :ceigen-lite))
    (when (or (and (not (fboundp s))
                   (not (boundp s)))
              (and (fboundp s)
                   (alexandria:starts-with-subseq "__" (symbol-name s))))
      (unintern s :ceigen-lite))))

(in-package :ceigen-lite)

(cl:let ((shared-library-pathname (cl:merge-pathnames #P"ceigen_lite/libceigen_lite.so" *src-dir*)))
  (cl:unless (cl:probe-file shared-library-pathname)
    (uiop:run-program
     (uiop:strcat "cd "
                  "'"
                  (cl:namestring (cl:merge-pathnames #P"ceigen_lite/" *src-dir*))
                  "'"
                  " && bash '"
                  (cl:namestring (cl:merge-pathnames #P"ceigen_lite/make.sh" *src-dir*))
                  "'")
     :error-output cl:t))
  (cffi:load-foreign-library shared-library-pathname))
