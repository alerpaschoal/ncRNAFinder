
#!/bin/bash

NORMAL="\\033[0;39m"
RED="\\033[0;31m"
GREEN="\\033[0;32m"
YELLOW="\\033[0;33m"
BLUE="\\033[0;34m"
TEST=$(pwd)
SAVE=$(pwd) 

die() {
    echo "$RED""Exit - ""$*""$NORMAL" 1>&2
    exit 1
}

echo "$YELLOW""Make sure your internet connection works for your shell prompt under current user's privilege ...""$NORMAL";
echo "$BLUE""Starting ncRNA-finder installation ...""$NORMAL";

#check for make
which make > /dev/null;
if [ $? != "0" ]; then
	echo "$RED""Can not proceed without make, please install and re-run ""$NORMAL"
	exit 1;
fi

#check for g++
which g++ > /dev/null;
if [ $? != "0" ]; then
	echo "$RED""Can not proceed without g++, please install and re-run""$NORMAL"
	exit 1;
fi

# check for gzip 
which gzip > /dev/null;
if [ $? != "0" ]; then
    echo "$RED""Can not proceed without gzip, please install and re-run""$NORMAL"
    exit 1;
fi

#check OS (Unix/Linux or Mac)
os=`uname`;

# get the right download program
if [ "$os" = "Darwin" ]; then
	# use curl as the download program 
	get="curl -L -o"
else
	# use wget as the download program
	get="wget --no-check-certificate -nc -O"
fi

if [ -d ./tmp_ncRNA-finder ]; then
    rm -rf ./tmp_ncRNA-finder
fi
mkdir ./tmp_ncRNA-finder
cd ./tmp_ncRNA-finder

################ Install dependencies  ###################

PREFIX_BIN=/usr/local/bin

if [ ! -w $PREFIX_BIN ]; then
    die "Cannot write to directory $PREFIX_BIN. Maybe missing super-user (root) permissions to write there.";
fi 

############### Install scritps ########################


chmod 777 $TEST/ncRNAFinder.py
chmod 777 $TEST/UpdateFasta.py
chmod 777 $TEST/Pipeline.py
chmod 777 $TEST/RNAcentral.py
chmod 777 $TEST/Rfam.py
echo "export PATH=\$PATH:$TEST" > /etc/profile.d/install_ncRNAFinder.sh
chmod 755 /etc/profile.d/install_ncRNAFinder.sh
CURRENT_PATH=$(cat /etc/environment | grep PATH | cut -d'=' -f2-)
echo "PATH=\"$CURRENT_PATH:$TEST\"" | sudo tee -a /etc/environment > /dev/null


which $TEST/ncRNAFinder.py > /dev/null;
if [ $? = "0" ]; then
	echo "$BLUE""ncRNA-finder appears to be installed successfully""$NORMAL"
else
	echo "$RED""ncRNA-finder NOT installed successfully""$NORMAL"; exit 1;
fi

################ Python Libraries ###################
echo "$BLUE""Installing Python Libraries ...""$NORMAL";

wasInstalled=0;
which pip > /dev/null;
if [ $? != "0" ]; then
	echo -n "$RED""Pip not installed. \n""$NORMAL"
	echo -n "$BLUE""Installing pip ... \n""$NORMAL"
	apt install python3-pip
	wasInstalled=1;
fi

python3 -m pip install -r ../requirements.txt



################ Infernal ###################

wasInstalled=0;
which cmscan > /dev/null;
if [ $? = "0" ]; then
	echo -n "$BLUE""Infernal appears to be already installed. \n""$NORMAL"
	wasInstalled=1;
fi

echo -n "Would you like to install/re-install Infernal? (y/n) [n] : "
read ans
if [ XX${ans} = XXy ]; then
	$get infernal-1.1.tar.gz http://eddylab.org/infernal/infernal-1.1.5-linux-intel-gcc.tar.gz
	tar xvzf infernal-1.1.tar.gz
	cp infernal-1.1.5-linux-intel-gcc/binaries/* $PREFIX_BIN
	wasInstalled=0;
fi

#some Infernal tests
if [ $wasInstalled = 0 ]; then
    which cmscan > /dev/null;
    if [ $? = "0" ]; then
	echo "$BLUE""Infernal appears to be installed successfully""$NORMAL"
    else
	echo "$RED""Infernal NOT installed successfully""$NORMAL"; exit 1;
    fi
fi

################ Rfam ###################

echo -n "Would you like to download/re-download Rfam database? (y/n) [y] : "
read ans
if [ XX${ans} != XXn ]; then
	echo "Where should Rfam be downloaded ? [$TEST]\nWrite new directory or press Enter to continue... "
	read answ
	answ=${answ:-$TEST}
	TEST=$answ
	if [ ! -d $TEST ]; then
		echo "Directory $TEST does not exist!"
    		echo -n "Do you want to create $TEST folder ? (y/n) [n] : "
    		read ans
    		if [ XX${ans} = XXy ]; then
        		mkdir $TEST || die "Cannot create  $TEST folder. Maybe missing super-user (root) permissions"
    		else
        		die "Must specify a directory to download Rfam.cm!"
    		fi
	fi
	if [ ! -w $TEST ]; then
    		die "Cannot write to directory $TEST. Maybe missing super-user (root) permissions to write there.";
	fi
	$get Rfam.clanin https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin 
	$get Rfam.cm.gz https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
	$get family.txt.gz https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/family.txt.gz
	gzip -d Rfam.cm.gz
	gzip -d family.txt.gz
	cmpress -F Rfam.cm
	cp Rfam.cm* $TEST
	cp Rfam.clanin $TEST
	cp family.txt $TEST
	$SAVE/Rfam.py
	
	echo "$BLUE""Rfam.cm appears to be downloaded successfully""$NORMAL"	
fi


################ BLASTn ################

wasInstalled=0;
which blastn > /dev/null;
if [ $? = "0" ]; then
	echo -n "$BLUE""BLAST appears to be already installed. \n""$NORMAL"
	wasInstalled=1;
fi

echo -n "Would you like to install/re-install BLAST? (y/n) [n] : "
read ans
if [ XX${ans} = XXy ]; then
	$get blast-2.15.0+-x64-linux.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.15.0+-x64-linux.tar.gz
	tar xvzf blast-2.15.0+-x64-linux.tar.gz
	cp ncbi-blast-2.15.0+/bin/* $PREFIX_BIN
	wasInstalled=0;
fi

#some Infernal tests
if [ $wasInstalled = 0 ]; then
    which blastn > /dev/null;
    if [ $? = "0" ]; then
	echo "$BLUE""BLAST appears to be installed successfully""$NORMAL"
    else
	echo "$RED""BLAST NOT installed successfully""$NORMAL"; exit 1;
    fi
fi




############### RNAcentral ###################
rnacentralBack()
{
$get rnacentral_active.fasta.gz https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/rnacentral_active.fasta.gz
$get id_mapping.tsv.gz https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/id_mapping/id_mapping.tsv.gz
gzip -d rnacentral_active.fasta.gz
gzip -d id_mapping.tsv.gz
cp rnacentral_active.fasta $TEST
cp id_mapping.tsv $TEST
$SAVE/RNAcentral.py
echo "$BLUE""rnacentral_active.fasta appears to be downloaded successfully""$NORMAL"

# clean up
chmod 777 $SAVE/*
cd ..
rm -rf tmp_ncRNA-finder

echo;
echo "$GREEN""Installation completed! Run ncRNAFinder.py in the terminal to check up the help guide.""$NORMAL"
}

echo -n "Would you like to download/re-download rnacentral database?""$RED""(Time-consuming step ~7 hours)""$NORMAL"" (y/n) [y] : "
read ans
if [ XX${ans} != XXn ]; then
	echo "Where should rnacentral_active.fasta be downloaded ? [$TEST]\nWrite new directory or press Enter to continue... "
	read answ
	answ=${answ:-$TEST}
	TEST=$answ
	if [ ! -d $TEST ]; then
		echo "Directory $TEST does not exist!"
    		echo -n "Do you want to create $TEST folder ? (y/n) [n] : "
    		read ans
    		if [ XX${ans} = XXy ]; then
        		mkdir $TEST || die "Cannot create  $TEST folder. Maybe missing super-user (root) permissions"
    		else
        		die "Must specify a directory to download RNAcentral.fasta!"
    		fi
	fi
	if [ ! -w $TEST ]; then
    		die "Cannot write to directory $TEST. Maybe missing super-user (root) permissions to write there.";
	fi
	
	echo -n "Would you like to run in background ? (y/n) [y] : " 
	read resp
	if [ XX${resp} = XXy ]; then
		rnacentralBack &
	else
		$get rnacentral_active.fasta.gz https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/rnacentral_active.fasta.gz
		$get id_mapping.tsv.gz https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/id_mapping/id_mapping.tsv.gz
		
		gzip -d rnacentral_active.fasta.gz
		gzip -d id_mapping.tsv.gz
		cp rnacentral_active.fasta $TEST
		cp id_mapping.tsv $TEST
		$SAVE/RNAcentral.py	

		echo "$BLUE""rnacentral_active.fasta appears to be downloaded successfully""$NORMAL"	
	



################ End of the installation ###################
		# clean up
		chmod 777 $SAVE/*
		cd ..
		rm -rf tmp_ncRNA-finder

		echo;
		echo "$GREEN""Installation completed! Run ncRNAFinder.py in the terminal to check up the help guide.""$NORMAL"
	fi
fi
