======
 AHF2
======
to include the contributions from a cube to its patch the following routines required modification:

+ patch_create()
+ patch_init()
+ patch_clear()
  => this is being taken care of by init_patch_physics() now

+ patch_add_psubcube(): add the actual subcube to patch
+ patch_include_adjacent_subcubes(): recursively search for neighbouring subcubes to be added to the same patch
  => this is being taken care of by add_patch_physics() now
  => note: patch_add_psubcube() is being called from inside patch_include_adjacent_subcubes()


TODO:
=====
- psc - ipart needs to point to memory address (and not position in array)
- properly normalize the patch->centre values
- generate patch-tree