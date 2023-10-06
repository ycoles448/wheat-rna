PROJ := "gosplit"
BIN := $(PROJ)

CC := gcc
LD := gcc

CFLAGS := $(CFLAGS) -std=c90 -Wall
LDFLAGS := $(LDFLAGS)

SRCDIR := .
OUTDIR := .

SRC := $(wildcard $(SRCDIR)/*.c)
OBJ := $(patsubst %.c,%.o,$(SRC))

.PHONY: clean default

$(OBJ): $(SRC)

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $<

$(BIN): $(OBJ)
	$(LD) $(LDFLAGS) -o $(BIN) $(OBJ)

debug: export CFLAGS += -g -DBUILD=debug
debug: $(BIN)

clean:
	rm -rf $(OBJ) $(BIN)
