package = "earth-magnetic-field"
version = "dev-1"
source = {
	url = "git+https://github.com/thenumbernine/earth-magnetic-field-lua.git"
}
description = {
	summary = "Visualization of the WMM 2020 data.",
	detailed = "Visualization of the WMM 2020 data.",
	homepage = "https://github.com/thenumbernine/earth-magnetic-field-lua",
	license = "MIT"
}
dependencies = {
	"lua ~> 5.1"
}
build = {
	type = "builtin",
	modules = {
		["earth-magnetic-field.run"] = "run.lua"
		["earth-magnetic-field.reduce"] = "reduce.lua"
	}
}
