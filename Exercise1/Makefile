CC = gcc

CFLAGS = -Wall -O3 -march=native
INCDIR = inc
SRCDIR = src
OBJDIR = obj

main.x: $(OBJDIR)/main.o $(OBJDIR)/pgm_utils.o
	$(CC) $(CFLAGS) $^ -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c -I$(INCDIR) $^ -o $@

clean:
	rm *.x
	rm $(OBJDIR)/*.o