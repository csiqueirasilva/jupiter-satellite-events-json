var fs = require('fs')

var regex = /^sat (\d): (\w{3}) (start|end)[ ]*: (.*)$/gm;

function parseFile (file) {
	
	var data = fs.readFileSync("data/" + file, 'utf8');

	data = data.replace(/[^]*Final results:/gm, '');
	
	var ret = {};

	var match = regex.exec(data);

	while(match !== null) {
		var o = {};
		o.sat = match[1];
		o.type = match[2];
		o.event = match[3];
		var date = new Date(match[4] + " GMT");
		o.date = date.getTime();
		
		var idx = date.toISOString().slice(0,10).replace(/-/g,"");
		
		if(!ret[idx]) {
			ret[idx] = [];
		}
		
		ret[idx].push(o);
		
		match = regex.exec(data);
	}
	
	var outputFile = file.replace("TXT", "json");
	
	fs.writeFileSync("export/" + outputFile, JSON.stringify(ret));

	console.log("Created: " + outputFile);
	
}

fs.readdir("data/.", function(err, data) {
	
	if(err) {
		return console.log(err);
	}
	
	for(var i = 0; i < data.length; i++) {
		if(data[i].match(/^.*\.TXT$/gm) !== null) {
			console.log('reading file ' + data[i]);
			parseFile(data[i]);
		}
	}
	
});