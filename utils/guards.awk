#! /usr/bin/gawk -f

BEGINFILE {
    first = 1;
    replace = 0;
}

/ifndef/ {
    if (first) {
	"realpath " FILENAME | getline fname;
	guard = toupper(substr(fname, index(fname, "mln/")));
	gsub(/\W/, "_", guard);
	oldguard = substr($0, index($0, "ifndef") + 7)
	if (oldguard != guard)
	{
	    print "Should I replace \n" "\t" oldguard "\nby\n" "\t" guard " (y/n)?";
	    while (1) {
		getline answer < "-";
		if (answer == "y" ) {
		    replace = 1; break;
		} else if (answer == "n") {
		    break;
		}
		print "Please answer y/n !"
	    }
	}
    }
}

replace && $0 ~ oldguard {
    sub(oldguard, guard);
}

{
    print $0;
}