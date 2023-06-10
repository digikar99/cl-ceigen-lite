(cl:in-package :cl)

(do-symbols (s (find-package :ceigen-lite))
  (when (and (fboundp s)
             (not (member s '(ceigen-lite::ceigen-lite-c-to-lisp
                              ceigen-lite::starts-with-ceigen-lite))))
    (fmakunbound s)
    (proclaim `(inline ,s))))
