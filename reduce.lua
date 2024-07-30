local ffi = require 'ffi'
local assertindex = require 'ext.assert'.index
local class = require 'ext.class'
local vector = require 'ffi.cpp.vector-lua'
local template = require 'template'
local gl = require 'gl'
local GLPingPong = require 'gl.pingpong'
local GLGeometry = require 'gl.geometry'
local GLProgram = require 'gl.program'
local GLTex2D = require 'gl.tex2d'

local GLReduce = class()

-- smallest width or height before switching to CPU reduce
GLReduce.minSize = 4

--[[
args:
	tex = texture to reduce.  only required if pingpong isn't provided.
	geom = (optional) geometry to use
	vertexShader (optional) vertex shader to use
	pingpong = pingpong to use (optional) default = make new one
	op = operation to reduce
--]]
function GLReduce:init(args)
	self.tex = args.tex
	self.geometry = args.geometry or GLGeometry{
		mode = gl.GL_TRIANGLE_STRIP,
		vertexes = {
			data = {
				0, 0,
				1, 0,
				0, 1,
				1, 1,
			},
			dim = 2,
		},
	}

	self.vertexShader = args.vertexShader or GLProgram.VertexShader{
		version = 'latest',
		precision = 'best',
		code = [[
layout(location=0) in vec2 vertex;
out vec2 texcoordv;
void main() {
	texcoordv = vertex;
	gl_Position = vec4(vertex * 2. - 1., 0., 1.);
}
]],
	}

	self.pingpong = args.pingpong or self:makePingPong{
		tex = self.tex,
		fbo = args.fbo,
	}

	self.gpuop = assertindex(args, 'gpuop')
	self.cpuop = assertindex(args, 'cpuop')
	self.program = GLProgram{
		version = 'latest',
		precision = 'best',
		shaders = {self.vertexShader},
		fragmentCode = template([[
in vec2 texcoordv;
out vec4 fragColor;
uniform sampler2D tex;
uniform vec2 offset;
uniform vec2 viewportSize;
uniform vec2 texSize;
void main() {
	// frag coord proportional to the viewport size, so ...
	// `* viewportSize` brings us to texel integers plus 1/2 offset
	// `/ texSize` puts us in unit texel coords
	vec2 texcoord = texcoordv * viewportSize / texSize;

	fragColor = <?=op(
		op(
			'texture(tex, texcoord)',
			'texture(tex, texcoord + vec2(offset.x, 0.))'
		),
		op(
			'texture(tex, texcoord + vec2(0., offset.y))',
			'texture(tex, texcoord + offset)'
		)
	)?>;
}
]],
		{
			op = self.gpuop,
		}),
	}

	local GLSceneObject = require 'gl.sceneobject'
	self.sceneobj = GLSceneObject{
		program = self.program,
		geometry = self.geometry,
	}

	local tex = self.tex or self.pingpong:cur()
	local ctype = assertindex(require 'gl.types'.ctypeForGLType, tex.type)
	local channels = assertindex(require 'gl.tex'.channelsForFormat, tex.format)	-- TODO move this table in gl.types?
	if channels == 1 then
	elseif channels == 2
	or channels == 3
	or channels == 4
	then
		local suffix = assertindex(require 'vec-ffi.suffix', ctype)
		ctype = require('vec-ffi.vec'..channels..suffix).name
	else
		error("idk what vector type to use for this many channels")
	end

	self.resultBuf = vector(ctype)
end

--[[ static helper
args:
	tex = tex to base off of
	fbo = fbo to use (optional)
--]]
function GLReduce:makePingPong(args)
	local tex = assertindex(args, 'tex')
	return GLPingPong{
		fbo = args.fbo,
		width = tex.width,
		height = tex.height,
		internalFormat = tex.internalFormat,
		format = tex.format,
		type = tex.type,
		minFilter = gl.GL_NEAREST,
		magFilter = gl.GL_NEAREST,
	}
end

function GLReduce:__call(tex)
--print('GLReduce:__call', self.gpuop('a', 'b'))
	tex = tex or self.tex
	assert(tex, "need a source tex")

	local sceneobj = self.sceneobj	-- will hold the vao or the attrs
	local pingpong = self.pingpong
	local fbo = pingpong.fbo
	fbo:bind()
	--[[ option #2) copy to the other texture
	fbo:setColorAttachmentTex2D(pingpong:cur().id)
	local res,err = fbo.check()
	if not res then print(err..'\n'..debug.traceback()) end
	srcTex:bind()
	--]]

	local program = self.program
	program:use()

-- TODO track this, maybe as a vector, maybe resize it, idk
	local resultBuf
	local srcTex = tex
	local srcW = srcTex.width
	local srcH = srcTex.height
	local dstW = math.ceil(srcW / 2)
	local dstH = math.ceil(srcH / 2)
	while true do 	--while srcW > self.minSize and srcH > self.minSize do
		local nextdone = srcW <= self.minSize and srcH <= self.minSize
--print('srcsize', srcW, srcH, 'dstsize', dstW, dstH, 'texsize', srcTex.width, srcTex.height)

		-- [[ option #1) change bound texture
		fbo:setColorAttachmentTex2D(pingpong:cur().id)
		local res,err = fbo.check()
		if not res then print(err..'\n'..debug.traceback()) end
		srcTex:bind()
		--]]

		gl.glViewport(0, 0, dstW, dstH)
--print('reduce offset', srcW - dstW, srcH - dstH)
		gl.glUniform2f(program.uniforms.offset.loc, (srcW - dstW) / srcTex.width, (srcH - dstH) / srcTex.height)
		gl.glUniform2f(program.uniforms.viewportSize.loc, dstW, dstH)
		gl.glUniform2f(program.uniforms.texSize.loc, srcTex.width, srcTex.height)

		if sceneobj.vao then
			sceneobj.vao:bind()
		else
			for _,attr in pairs(sceneobj.attrs) do
				attr:enableAndSet()
			end
		end
		sceneobj.geometry:draw()
		if sceneobj.vao then
			sceneobj.vao:unbind()
		else
			for _,attr in pairs(sceneobj.attrs) do
				attr:disable()
			end
		end

		-- if it's the last one ...
		-- TODO while fbo is bound
		-- but TODO dont do the reduce if the intiial size is too small
		if nextdone then
			self.resultBuf:resize(dstW * dstH)
			gl.glReadPixels(0, 0, dstW, dstH, tex.format, tex.type, self.resultBuf.v)
		end

		if nextdone then
			for i=1,dstW * dstH-1 do
				self.cpuop(self.resultBuf.v[0], self.resultBuf.v[i], self.resultBuf.v[0])
			end
			break
		end
		-- [[ option #1) change bound texture
		pingpong:swap()
		--]]
		--[[ option #2) copy to the other texture
		if it's the first pass then pingpong:prev():bind() end	 ... or just do a copy up front, and bind once ...
		gl.glCopyTexSubImage2D(self.fieldLinePosTex.target, 0, 0, 0, i, 0, 1, self.fieldLineVtxsTex.height)
		--]]
		srcTex = pingpong:prev()
		srcW = dstW
		srcH = dstH
		dstW = math.ceil(srcW / 2)
		dstH = math.ceil(srcH / 2)
	end

	program:useNone()
	fbo:unbind()
	GLTex2D:unbind()

	return self.resultBuf.v[0]
end

return GLReduce
