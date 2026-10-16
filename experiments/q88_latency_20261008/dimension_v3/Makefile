LAZER_ROOT ?= ../../../../lazer
NATIVE := $(LAZER_ROOT)/src/labrados
SOURCES := aesctr comkey constraints cpucycles data fips202 gaussian jlproj labrador_core labrador_tail labrador lnp ntt pack poly polz proofsystem randombytes rejection labrados_python timing
FLAGS := -std=gnu2x -O3 -DNDEBUG -fwrapv -ffast-math -fopenmp -flto=auto -fno-semantic-interposition -fno-plt -mavx512f -mavx512bw -mavx512dq -mavx512vl -mavx512vbmi -mavx512vbmi2 -mavx512vnni -mvaes -maes -mpclmul -fvisibility=hidden
LIBRARY ?= liblabrador38.so
$(LIBRARY): dachshund.c labradoodle.c polx.c scratch_pool.h $(addprefix $(NATIVE)/,$(addsuffix .c,$(SOURCES)))
	$(CC) $(FLAGS) -DLOGQ=38 -I$(NATIVE) -I. -shared -fPIC -o $@ dachshund.c labradoodle.c polx.c $(addprefix $(NATIVE)/,$(addsuffix .c,$(SOURCES))) -lmvec -lm
