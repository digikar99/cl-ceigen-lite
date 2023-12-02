(cl:in-package :cl)

(eval-when (:compile-toplevel :load-toplevel :execute)
  (do-symbols (s (find-package :ceigen-lite))
    (when (or (and (not (fboundp s))
                   (not (boundp s)))
              (and (fboundp s)
                   (alexandria:starts-with-subseq "__" (symbol-name s))))
      (unintern s :ceigen-lite))))

(in-package :ceigen-lite)

(cl:let ((shared-library-pathname
           (cl:merge-pathnames (cl:pathname
                                (cl:format cl:nil "ceigen_lite/libceigen_lite-~A-~A"
                                           (cl:substitute #\_ #\- (cl:string-downcase (cl:machine-type)))
                                           (cl:string-downcase
                                            (cl:cond ((cl:member :linux cl:*features*)
                                                      "linux.so")
                                                     ((cl:member :windows cl:*features*)
                                                      "windows.dll")
                                                     (cl:t
                                                      (cl:error "Unhandled operating system in cl-ceigen-lite"))))))
                               *src-dir*)))
  (cl:unless (cl:probe-file shared-library-pathname)
    (uiop:run-program
     (uiop:strcat "cd "
                  "'"
                  (cl:namestring (cl:merge-pathnames #P"ceigen_lite/" *src-dir*))
                  "'"
                  " && bash '"
                  (cl:namestring (cl:merge-pathnames #P"ceigen_lite/make.sh" *src-dir*))
                  "'")
     :output cl:t
     :error-output cl:t))
  (cffi:load-foreign-library shared-library-pathname))
