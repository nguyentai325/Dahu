#! /usr/bin/gawk -f

BEGINFILE {
    inbox = 0;
}



/\/\*\*\*\*+\// && !inbox {
    inbox=1;
}


/\/\*.+\*+\// {
    if (inbox) {
	match($0, /\*+\s+(.+\w+)\*/, arr);
	label=arr[1];
    };
}

/\/\*\*\*\*+\// && !inbox {
    inbox=1;
}



