def materials():
	return {
		'ceramic': {
			'youngModulus': 7.4e10,
			'poissonRatio': 0.19,
			'density': 2300.0,
			'alpha': 6.0,
			'beta': 1e-7,
			'gamma': 3e-2,
			'restitution_coeff': 0.4,
			'friction_coeff': 0.2
		},
		'polystyrene': {
			'youngModulus': 3.5e9,
			'poissonRatio': 0.34,
			'density': 1050.0,
			'alpha': 30.0,
			'beta': 8e-7,
			'gamma': 4e-4,
			'restitution_coeff': 0.4,
			'friction_coeff': 0.2
		},
		'steel': {
			'youngModulus': 2e11,
			'poissonRatio': 0.29,
			'density': 7850.0,
			'alpha': 5.0,
			'beta': 3e-8,
			'gamma': 9e-3,
			'restitution_coeff': 0.4,
			'friction_coeff': 0.2
		},
		'mdf': {
			'youngModulus': 4e9,
			'poissonRatio': 0.32,
			'density': 615.0,
			'alpha': 35.0,
			'beta': 5e-6,
			'gamma': 9e-3,
			'restitution_coeff': 0.4,
			'friction_coeff': 0.2
		},
		'wood': {
			'youngModulus': 1.1e10,
			'poissonRatio': 0.25,
			'density': 750.0,
			'alpha': 60.0,
			'beta': 2e-6,
			'gamma': 5e-4,
			'restitution_coeff': 0.4,
			'friction_coeff': 0.2
		},		
	}