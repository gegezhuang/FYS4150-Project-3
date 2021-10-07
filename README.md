# FYS4150 Project 3

# Usage

## During development
To debug, use `make debug`. This will compile main.cpp` with some extra flags, and open gdb to see where things went wrong.

If you want to test all the tests in `test.cpp`, run `make test`.

If you want to compile the pdf, use `make project.pdf`, and if you want to force everything to remake, just use `make`.

When you have made something that needs to run before compiling `project.pdf`, add it to the `ALL_TARGET_FILES` variable in the makefile, and then make a rule to make the target file at the bottom of the makefile.

## After completion
```bash
make
```
