#== makefile for the layer ==

all:
	@echo ================================================================
	@echo ================== Building CLC layer ==========================
	@echo ================================================================
	make -C bpl  -f bpl.mak bits=$(bits) compiler=$(compiler)
	make -C ipl  -f ipl.mak bits=$(bits) compiler=$(compiler)

clean:
	@echo ================================================================
	@echo ================== Cleaning CLC layer ==========================
	@echo ================================================================
	make -C bpl  -f bpl.mak clean bits=$(bits) compiler=$(compiler)
	make -C ipl  -f ipl.mak clean bits=$(bits) compiler=$(compiler)
