(cl:in-package :cl)

(eval-when (:compile-toplevel :load-toplevel :execute)
  (do-symbols (s (find-package :ceigen-lite))
    (when (or (and (not (fboundp s))
                   (not (boundp s)))
              (and (fboundp s)
                   (alexandria:starts-with-subseq "__" (symbol-name s))))
      (unintern s :ceigen-lite))))

(in-package :ceigen-lite)

(cl:defparameter *ceigen-lite-release-tag* "v3.4-2026.07"
  "Update this when a new ceigen_lite release is cut.")

(cl:let* ((arch (cl:cond ((cl:member :x86-64 cl:*features*) "x86_64")
                         ((cl:member :arm64 cl:*features*) "arm64")
                         (cl:t (cl:error "cl-ceigen-lite: unrecognized CPU feature in ~S; ~
                                           don't know which ceigen_lite release asset to fetch."
                                         cl:*features*))))
          (os-ext (cl:cond ((cl:member :linux cl:*features*) "linux.so")
                           ((cl:member :darwin cl:*features*) "macos.dylib")
                           ((cl:member :windows cl:*features*) "windows.dll")
                           (cl:t (cl:error "Unhandled operating system in cl-ceigen-lite"))))
          (filename (cl:format cl:nil "libceigen_lite-~A-~A" arch os-ext))
          (shared-library-pathname (cl:merge-pathnames (cl:pathname (cl:format cl:nil "ceigen_lite/~A" filename))
                                                       *src-dir*))
          (download-url (cl:format cl:nil "https://github.com/digikar99/ceigen_lite/releases/download/~A/~A"
                                   *ceigen-lite-release-tag* filename)))
  (cl:unless (cl:probe-file shared-library-pathname)
    (cl:ensure-directories-exist shared-library-pathname)
    (cl:format cl:t "~&cl-ceigen-lite: downloading ~A~%" download-url)
    (uiop:run-program (cl:list "curl" "-L" "--fail" "-o" (cl:namestring shared-library-pathname) download-url)
                      :output cl:t
                      :error-output cl:t))
  (cffi:load-foreign-library shared-library-pathname))
