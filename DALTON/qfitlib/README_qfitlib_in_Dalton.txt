14. Febuary 2025 Hans Jørgen Aagaard Jensen

After consulting with the author of qfitlib, Casper Steinmann,
we agreed that I could include qfitlib directly in Dalton instead of 
as an external module, in order to make the inclusion more robust.

The initial version of today is a copy of git hash e336494 from
  https://github.com/cstein/qfitlib.git

Note that the PRG_DALTON must be defined.
----
Note that CMakeList.txt here is not used,
what is compiled is determined in the ../../cmake/SourcesDALTON.cmake file.
