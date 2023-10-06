PROJ := star2table
BIN := $(PROJ)

CC := gcc
LD := gcc

CFLAGS := $(CFLAGS) -std=c99 -Wall -Werror
LDFLAGS := $(LDFLAGS)

SRCDIR := src
SRC := $(wildcard $(SRCDIR)/*.c)
OBJ := $(patsubst %.c,%.o,$(SRC))

.PHONY: default clean

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $<

$(OBJ): $(SRC)

$(BIN): $(OBJ)
	$(LD) $(LDFLAGS) -o $(BIN) $(OBJ)

release: export CFLAGS += -DRELEASE
release: $(BIN)

debug: export CFLAGS += -DDEBUG -g
debug: $(BIN)

default: debug

clean:
	rm -rf $(BIN) $(OBJ)
