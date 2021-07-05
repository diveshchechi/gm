//%******[Time increase]*************************************************
	if (counting < 1003) {
		datsave(counting,lloyd);
	}
	else if (counting % 5 == 0) {
		datsave(counting,lloyd);
	}
	if ((counting-1000)%250 == 0) {
		// calculate new position of seeds
		findNewPosition();
		lloyd++;
		//counting = 1000;
		counting++;
		datsave(counting, lloyd);
		if (lloyd == 20) return 0;
		goto start;
		return 0;
	}
	counting++;
	goto start;
	return 0;
