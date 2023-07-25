#### SETUP FILE for you to prepare a CSV file with a Header and 5 columns, in the following order:

* Filename of Loupe browser DEG output file to be converted to SEGEX input format (filename without csv extension)
* Condition2 **Label for expression data to be used as the Numerator (<span style="color:red">Treatment</span>). Use text without spaces**.
* Condition1 **Label for expression data to be used as the Denominator (<span style="color:red">Control</span>). Use text without spaces**.
* Datatype identifier (**IntronicMonoExonic** or **Genebody**)
* ReverseColumns (0 or 1 where 0 means - do nothing and 1  means - reverse conditions; see instructions below)

DEG data exported from the Loupe browser (csv file) contains expression data for both conditions (<span style="color:red">Treatment, Control</span>). We need to manually specify which expression data is Treatment and which data is Control. We use the convention of Condition2 (numerator) to designate Treatment and Condition1 (denominator) to designate the Control.