# Legacy reference implementations

The source workspace contains two legacy packages used for provenance and regression validation:

- `cellulose-builder-2013/`: Cellulose-Builder by Thiago C. F. Gomes, licensed under GPLv3. Cite Gomes and Skaf, *Journal of Computational Chemistry* 33, 1338-1346 (2012), DOI: 10.1002/jcc.22959.
- `doglycans-2017/`: doGlycans by Reinis Danne and collaborators, licensed under GPLv3. Cite Danne et al., *Journal of Chemical Information and Modeling* 57, 2401-2406 (2017), DOI: 10.1021/acs.jcim.7b00237.
- `source-archives/`: unmodified source and supporting-information archives.
- `publications/`: publications retained for scientific reference.

These directories are retained unchanged as reference material. The new package does not import them at runtime. Any future compatibility patch must be isolated, documented, and retain the original copyright and license notices.
