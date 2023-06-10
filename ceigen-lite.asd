(defsystem "ceigen-lite"
  :pathname ""
  :author "Shubhamkar B. Ayare (shubhamayare@yahoo.co.in)"
  :description "A Common Lisp wrapper around CEIGEN-LITE - which is itself a C wrapper around the C++ Eigen library."
  :license "MIT"
  :depends-on ("uiop"
               "cffi"
               "cl-autowrap")
  :serial t
  :components ((:module "specs"
                :components ((:static-file "ceigen_lite.aarch64-pc-linux-gnu.spec")
                             (:static-file "ceigen_lite.aarch64-unknown-linux-android.spec")
                             (:static-file "ceigen_lite.arm-pc-linux-gnu.spec")
                             (:static-file "ceigen_lite.arm-unknown-linux-androideabi.spec")
                             (:static-file "ceigen_lite.i386-unknown-freebsd.spec")
                             (:static-file "ceigen_lite.i386-unknown-openbsd.spec")
                             (:static-file "ceigen_lite.i686-apple-darwin9.spec")
                             (:static-file "ceigen_lite.i686-pc-linux-gnu.spec")
                             (:static-file "ceigen_lite.i686-pc-windows-msvc.spec")
                             (:static-file "ceigen_lite.i686-unknown-linux-android.spec")
                             (:static-file "ceigen_lite.x86_64-apple-darwin9.spec")
                             (:static-file "ceigen_lite.x86_64-pc-linux-gnu.spec")
                             (:static-file "ceigen_lite.x86_64-pc-windows-msvc.spec")
                             (:static-file "ceigen_lite.x86_64-unknown-freebsd.spec")
                             (:static-file "ceigen_lite.x86_64-unknown-linux-android.spec")
                             (:static-file "ceigen_lite.x86_64-unknown-openbsd.spec")))
               (:module "src"
                :serial t
                :components ((:file "package")
                             (:file "noninline")
                             ;; FIXME: Simplify after autowrap adds an inline option
                             (:file "inline-declarations")
                             (:file "inline")
                             (:file "shared-object")))))
