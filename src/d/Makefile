D_OPTIMIZE=1
D_USE_FREELISTS=1

MAJOR=1
MINOR=30
RELEASE=$(MAJOR).$(MINOR)


LIB_SRCS = arg.c parse.c scan.c symtab.c util.c dparse_tree.c
LIB_OBJS = $(LIB_SRCS:%.c=%.o)

LIBDPARSE = libdparse.a

$(LIBDPARSE): $(LIB_OBJS)
	ar crv $@ $^


