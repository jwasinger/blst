.text

.p2align	5
.Lone_256:
.quad	1,0,0,0

.globl	eucl_inverse_mod_256

.def	eucl_inverse_mod_256;
.type	32;
.endef
.p2align	5
eucl_inverse_mod_256:
.long	3573752639
	stp	x29,x30,[sp,#-16]!
	add	x29,sp,#0
	sub	sp,sp,#128

	adr	x4,.Lone_256
	cmp	x3,#0
	csel	x3,x3,x4,ne		// default x3 to 1

	ldp	x4,x5,[x1]
	ldp	x6,x7,[x1,#16]

	orr	x8,x4,x5
	orr	x9,x6,x7
	orr	x10,x8,x9
	cbz	x10,.Labort_256		// abort if |inp|==0

	add	x1, sp, #0

	ldp	x8,x9,[x3]
	ldp	x10,x11,[x3,#16]

	ldp	x12,x13,[x2]
	ldp	x14,x15,[x2,#16]

	add	x2, sp, #64

	stp	x4,x5,[x1,#0]	// copy |inp| to U
	stp	x6,x7,[x1,#16]

	stp	x8,x9,[x1,#32]	// copy |one| to X
	stp	x10,x11,[x1,#48]

	stp	x12,x13,[x2,#0]	// copy |mod| to V
	stp	x14,x15,[x2,#16]

	stp	xzr,xzr,[x2,#32]		// clear Y
	stp	xzr,xzr,[x2,#48]
	b	.Loop_inv_256

.p2align	4
.Loop_inv_256:
	####### *x2 is always odd at this point, no need to check...

	bl	__remove_powers_of_2_256

	ldp	x8,x9,[x2]
	ldp	x10,x11,[x2,#16]
	subs	x4,x4,x8		// U-V
	sbcs	x5,x5,x9
	sbcs	x6,x6,x10
	sbcs	x7,x7,x11
	b.hs	.Lu_greater_than_v_256

	eor	x2,x2,x1		// xchg	x2,x1
	mvn	x4,x4			// U-V => V-U
	eor	x1,x1,x2
	mvn	x5,x5
	eor	x2,x2,x1
	adds	x4,x4,#1
	mvn	x6,x6
	adcs	x5,x5,xzr
	mvn	x7,x7
	adcs	x6,x6,xzr
	adc	x7,x7,xzr

.Lu_greater_than_v_256:
	stp	x4,x5,[x1]
	orr	x3,x4,x5
	ldp	x4,x5,[x2,#32]
	orr	x3,x3,x6
	ldp	x8,x9,[x1,#32]
	orr	x3,x3,x7
	stp	x6,x7,[x1,#16]
	ldp	x6,x7,[x2,#48]
	subs	x8,x8,x4		// X-Y		# [alt. Y-X]
	ldp	x10,x11,[x1,#48]
	sbcs	x9,x9,x5
	sbcs	x10,x10,x6
	sbcs	x11,x11,x7
	sbc	x7,xzr,xzr			// borrow -> mask

	and	x4,x12,x7
	and	x5,x13,x7
	adds	x8,x8,x4		// reduce if X<Y # [alt. Y<X]
	and	x6,x14,x7
	adcs	x9,x9,x5
	and	x7,x15,x7
	adcs	x10,x10,x6
	stp	x8,x9,[x1,#32]
	adc	x11,x11,x7
	stp	x10,x11,[x1,#48]

	cbnz	x3,.Loop_inv_256		// U!=0?

	ldr	x30,[x29,#8]
	ldp	x4,x5,[x2,#32]	// return Y
	ldp	x6,x7,[x2,#48]
	mov	x10,#1

.Labort_256:
	stp	x4,x5,[x0]
	stp	x6,x7,[x0,#16]
	mov	x0,x10			// boolean return value

	add	sp,sp,#128
	ldr	x29,[sp],#16
.long	3573752767
	ret


.def	__remove_powers_of_2_256;
.type	32;
.endef
.p2align	4
__remove_powers_of_2_256:
	ldp	x4,x5,[x1]
	ldp	x6,x7,[x1,#16]

.Loop_of_2_256:
	rbit	x3,x4
	tbnz	x4,#0,.Loop_of_2_done_256

	clz	x3,x3
	cmp	x4,#0
	mov	x8,#63
	csel	x3,x3,x8,ne		// unlikely in real life

	neg	x11,x3
	lsrv	x4,x4,x3		// acc[0:3] >>= cnt
	lslv	x8,x5,x11
	orr	x4,x4,x8
	lsrv	x5,x5,x3
	lslv	x9,x6,x11
	orr	x5,x5,x9
	ldp	x8,x9,[x1,#32]
	lsrv	x6,x6,x3
	lslv	x10,x7,x11
	orr	x6,x6,x10
	ldp	x10,x11,[x1,#48]
	lsrv	x7,x7,x3

	stp	x4, x5,[x1]
	stp	x6, x7,[x1,#16]
	b	.Loop_div_by_2_256

.p2align	4
.Loop_div_by_2_256:
	sbfx	x7,x8,#0,#1
	sub	x3,x3,#1

	and	x4,x12,x7
	and	x5,x13,x7
	adds	x8,x8,x4
	and	x6,x14,x7
	adcs	x9,x9,x5
	and	x7,x15,x7
	adcs	x10,x10,x6
	extr	x8,x9,x8,#1	// acc[4:7] >>= 1
	adcs	x11,x11,x7
	extr	x9,x10,x9,#1
	adc	x7,xzr,xzr		// redundant if modulus is <256 bits...
	extr	x10,x11,x10,#1
	extr	x11,x7,x11,#1

	cbnz	x3,.Loop_div_by_2_256

	ldp	x4,x5,[x1] // reload X [mostly for 2nd caller]
	ldp	x6,x7,[x1,#16]

	stp	x8,x9,[x1,#32]
	stp	x10,x11,[x1,#48]

	tbz	x4,#0,.Loop_of_2_256// unlikely in real life

.Loop_of_2_done_256:
	ret

