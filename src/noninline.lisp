(cl:in-package :ceigen-lite)

(autowrap:c-include (cl:merge-pathnames #P"ceigen_lite/ceigen_lite.h" *src-dir*)
                    :spec-path
                    (cl:merge-pathnames #P"specs/" *src-dir*)
                    :c-to-lisp-function #'ceigen-lite-c-to-lisp
                    :include-definitions #'starts-with-ceigen-lite)
