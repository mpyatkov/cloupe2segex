#### CSV file with description (4 columns with header):

Contains cloupe csv filenames which should be converted to SEGEX format (tsv)

* Current filename without extension
* Condition 1 (**text without spaces**)
* Condition 2 (**text without spaces**)
* Identifier (**IntronicMonoExonic** or **Genebody**)
* ReverseColumns (**0** or **1** where 0 - do nothing and 1 - reverse conditions)

The DEG data exported from the Loupe browser contains information about both conditions (Control, Treatment) at the same time, and we have to manually specify what is Control and what is Treatment to make it easier to import into SEGEX. This service, like SEGEX, uses the Condition2 (numerator) and Condition1 (denominator) notation by default.
