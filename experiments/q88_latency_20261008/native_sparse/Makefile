LAZER_ROOT ?= ../../../../lazer
NATIVE := $(LAZER_ROOT)/src/labrados
SOURCES := aesctr comkey constraints cpucycles dachshund data fips202 gaussian jlproj labradoodle labrador_core labrador_tail labrador lnp ntt pack poly polz proofsystem randombytes rejection labrados_python timing
FLAGS := -std=gnu2x -O3 -DNDEBUG -fwrapv -ffast-math -fopenmp -flto=auto -fno-semantic-interposition -fno-plt -mavx512f -mavx512bw -mavx512dq -mavx512vl -mavx512vbmi -mavx512vbmi2 -mavx512vnni -mvaes -maes -mpclmul -fvisibility=hidden
liblabrador38.so: polx.c scratch_pool.h $(addprefix $(NATIVE)/,$(addsuffix .c,$(SOURCES)))
	$(CC) $(FLAGS) -DLOGQ=38 -I$(NATIVE) -I. -shared -fPIC -o $@ polx.c $(addprefix $(NATIVE)/,$(addsuffix .c,$(SOURCES))) -lmvec -lm
