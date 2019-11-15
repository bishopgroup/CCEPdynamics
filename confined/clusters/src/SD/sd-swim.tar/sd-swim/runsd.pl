#!/usr/local/bin/perl

	for ( $i = 1; $i < 3; $i++ )
	{

		open( seed, ">seed.in" ) or die "Couldn't open file";

		$myseed = time ^ $$ ^ unpack "%L*", `ps axww | gzip `; 

                print seed "$myseed"; 

                close( seed );

		system( '../sd/./sd.x' );

                if ( $i < 10 ) {
	
			system( "cp pos.out analysis/pos000$i.out" );
			system( "cp drift.out analysis/drift000$i.out" );

		} elsif ( $i < 100 ) {

			system( "cp pos.out analysis/pos00$i.out" );
			system( "cp drift.out analysis/drift00$i.out" );

		} elsif ( $i < 1000 ) {

			system( "cp pos.out analysis/pos0$i.out" );
			system( "cp drift.out analysis/drift0$i.out" );

		} else {

			system( "cp pos.out analysis/pos$i.out" );
			system( "cp drift.out analysis/drift$i.out" );
		
		}

	}


