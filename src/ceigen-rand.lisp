(defpackage :ceigen-rand
  (:use :cl)
  (:import-from #:ceigen-lite #:seed)
  (:export

   #:seed

   #:beta
   #:beta!
   #:cauchy
   #:cauchy!
   #:chi-squared
   #:chi-squared!
   #:exponential
   #:exponential!
   #:extreme-value
   #:extreme-value!
   #:fisher-f
   #:fisher-f!
   #:gamma
   #:gamma!
   #:lognormal
   #:lognormal!
   #:normal
   #:normal!
   #:uniform-real
   #:uniform-real!
   #:weibull
   #:weibull!

   #:bernoulli
   #:bernoulli!
   #:binomial
   #:binomial!
   #:geometric
   #:geometric!
   #:negative-binomial
   #:negative-binomial!
   #:poisson
   #:poisson!
   #:uniform-int
   #:uniform-int!))

(in-package :ceigen-rand)

(defmacro define-random-fn (name &rest params)
  "
NAME     - symbol naming the distribution, e.g. NORMAL, BETA, CHISQUARE.
PARAMS   - list of (param-name default-value), e.g. ((mean 0) (sd 1)).
           Order must match the CEIGEN-LITE C function's argument order.
Expects CEIGEN-LITE to export S<NAME> and D<NAME> with signature
  (length pointer &rest params) for single-float / double-float respectively."
  (let* ((default-type 'double-float)
         (symbol-name (symbol-name name))
         (bang-name  (alexandria:symbolicate name "!"))
         (param-syms (mapcar #'first params))
         (s-fn (find-symbol (format nil "S~A" symbol-name) :ceigen-lite))
         (d-fn (find-symbol (format nil "D~A" symbol-name) :ceigen-lite))
         (i32-fn (find-symbol (format nil "I32~A" symbol-name) :ceigen-lite))
         (type (gensym "TYPE")))
    (when i32-fn (setf default-type '(signed-byte 32)))
    (if (equal default-type '(signed-byte 32))
        (unless (and i32-fn (fboundp i32-fn))
          (warn "~A: CEIGEN-LITE:~A not found — ceigen-rand may not implement ~
             this distribution yet."
                name
                (format nil "I32~A" symbol-name)))
        (unless (and s-fn d-fn (fboundp s-fn) (fboundp d-fn))
          (warn "~A: CEIGEN-LITE:~A / :~A not found — ceigen-rand may not implement ~
             this distribution yet."
                name
                (format nil "S~A" symbol-name)
                (format nil "D~A" symbol-name))))
    `(progn
       (defun ,bang-name (vector &optional ,@params)
         (declare (type ,(if (equal default-type '(signed-byte 32))
                             `(simple-array (signed-byte 32) 1)
                             `(or (simple-array single-float 1)
                                  (simple-array double-float 1)))
                        vector)
                  #+sbcl (sb-ext:muffle-conditions sb-ext:compiler-note)
                  (optimize speed))
         (let* ((,type (array-element-type vector))
                ,@(mapcar (lambda (p)
                            `(,(first p) ,(cond ((integerp (second p))
                                                 (first p))
                                                ((equal default-type '(signed-byte 32))
                                                 `(coerce ,(first p) 'double-float))
                                                (t
                                                 `(coerce ,(first p) ,type)))))
                          params)
                (len (length vector)))
           (declare (ignorable ,type ,@param-syms))
           (cffi:with-pointer-to-vector-data (ptr vector)
             ,(if (equal default-type '(signed-byte 32))
                  `(etypecase vector
                     ((vector (signed-byte 32)) (,i32-fn len ptr ,@param-syms)))
                  `(etypecase vector
                     ((vector single-float) (,s-fn len ptr ,@param-syms))
                     ((vector double-float) (,d-fn len ptr ,@param-syms)))))
           vector))

       (defun ,name (length &optional ,@params (type ',default-type))
         (declare (optimize speed)
                  #+sbcl (sb-ext:muffle-conditions sb-ext:compiler-note)
                  (type fixnum length))
         (let ((array (make-array length :element-type type)))
           (,bang-name array ,@param-syms))))))

(define-random-fn beta        (a 1.0) (b 1.0))
(define-random-fn cauchy      (a 0.0) (b 1.0))
(define-random-fn chi-squared (df 1.0))
(define-random-fn exponential (lambda 1.0))
(define-random-fn extreme-value (loc 0.0) (scale 1.0))
(define-random-fn fisher-f    (m 1.0) (n 1.0))
(define-random-fn gamma       (alpha 1.0) (beta 1.0))
(define-random-fn lognormal  (mean 0.0) (stdev 1.0))
(define-random-fn normal      (mean 0.0) (sd 1.0))
(define-random-fn student-t   (ndof 1.0))
(define-random-fn uniform-real (min 0.0) (max 1.0))
(define-random-fn weibull     (a 1.0) (b 1.0))


(define-random-fn bernoulli         (p 0.5))
(define-random-fn binomial          (trials 1) (p 0.5))
(define-random-fn geometric         (p 0.5))
(define-random-fn negative-binomial (trials 1) (p 0.5))
(define-random-fn poisson           (mean 1.0))
(define-random-fn uniform-int       (min 0) (max 1))
