import React, { useEffect, useRef, useState } from "react"
import * as THREE from "three"

const TorusIrrationalLine = () => {
  const mountRef = useRef(null)
  const [slope, setSlope] = useState(Math.sqrt(2))
  const [maxWraps, setMaxWraps] = useState(100)
  const [isPlaying, setIsPlaying] = useState(true)
  const [speed, setSpeed] = useState(1)

  useEffect(() => {
    const scene = new THREE.Scene()
    scene.background = new THREE.Color(0x0a0a0a)

    const camera = new THREE.PerspectiveCamera(
      75,
      window.innerWidth / window.innerHeight,
      0.1,
      1000
    )
    camera.position.set(6, 4, 6)
    camera.lookAt(0, 0, 0)

    const renderer = new THREE.WebGLRenderer({ antialias: true })
    renderer.setSize(window.innerWidth, window.innerHeight)
    mountRef.current.appendChild(renderer.domElement)

    // Create torus with proper alignment
    const R = 2.5 // Major radius
    const r = 0.8 // Minor radius
    const torusGeometry = new THREE.TorusGeometry(R, r, 64, 128)
    const torusMaterial = new THREE.MeshStandardMaterial({
      color: 0x2196f3,
      transparent: true,
      opacity: 0.4,
      side: THREE.DoubleSide,
      roughness: 0.5,
      metalness: 0.3
    })
    const torus = new THREE.Mesh(torusGeometry, torusMaterial)
    torus.rotation.x = Math.PI / 2 // Rotate 90 degrees
    scene.add(torus)

    // Create wireframe overlay
    const wireframeGeometry = new THREE.TorusGeometry(R, r, 32, 64)
    const wireframeMaterial = new THREE.MeshBasicMaterial({
      color: 0x64b5f6,
      wireframe: true,
      transparent: true,
      opacity: 0.15
    })
    const wireframe = new THREE.Mesh(wireframeGeometry, wireframeMaterial)
    wireframe.rotation.x = Math.PI / 2 // Rotate 90 degrees
    scene.add(wireframe)

    // Create line geometry that will be updated
    const maxPoints = maxWraps * 200
    const positions = new Float32Array(maxPoints * 3)
    const lineGeometry = new THREE.BufferGeometry()
    lineGeometry.setAttribute(
      "position",
      new THREE.BufferAttribute(positions, 3)
    )
    lineGeometry.setDrawRange(0, 0)

    const lineMaterial = new THREE.LineBasicMaterial({
      color: 0xff5252,
      linewidth: 2
    })
    const line = new THREE.Line(lineGeometry, lineMaterial)
    line.rotation.x = Math.PI / 2 // Rotate 90 degrees
    scene.add(line)

    // Point at the current end of the line
    const pointGeometry = new THREE.SphereGeometry(0.08, 16, 16)
    const pointMaterial = new THREE.MeshBasicMaterial({ color: 0xffeb3b })
    const point = new THREE.Mesh(pointGeometry, pointMaterial)
    point.rotation.x = Math.PI / 2 // Rotate 90 degrees
    scene.add(point)

    // Lighting
    const ambientLight = new THREE.AmbientLight(0xffffff, 0.4)
    scene.add(ambientLight)

    const directionalLight1 = new THREE.DirectionalLight(0xffffff, 0.6)
    directionalLight1.position.set(5, 5, 5)
    scene.add(directionalLight1)

    const directionalLight2 = new THREE.DirectionalLight(0xffffff, 0.3)
    directionalLight2.position.set(-5, 3, -5)
    scene.add(directionalLight2)

    // Animation state
    let currentWrap = 0
    let autoRotate = 0

    // Precompute all points (adjusted for 90-degree rotated torus)
    const computePoint = (t, s) => {
      const theta = t
      const phi = s

      const x = (R + r * Math.cos(phi)) * Math.cos(theta)
      const y = (R + r * Math.cos(phi)) * Math.sin(theta)
      const z = r * Math.sin(phi)

      return { x, y, z }
    }

    // Animation loop
    let animationId
    const animate = () => {
      animationId = requestAnimationFrame(animate)

      if (isPlaying && currentWrap < maxWraps) {
        currentWrap += 0.02 * speed

        const numPoints = Math.floor(currentWrap * 200)
        const positionAttribute = lineGeometry.getAttribute("position")

        for (let i = 0; i <= numPoints && i < maxPoints; i++) {
          const t = (i / 200) * 2 * Math.PI
          const s = slope * t
          const pos = computePoint(t, s)

          positionAttribute.setXYZ(i, pos.x, pos.y, pos.z)
        }

        positionAttribute.needsUpdate = true
        lineGeometry.setDrawRange(0, Math.min(numPoints, maxPoints))

        // Update point position
        if (numPoints > 0) {
          const t = (numPoints / 200) * 2 * Math.PI
          const s = slope * t
          const pos = computePoint(t, s)
          point.position.set(pos.x, pos.y, pos.z)
        }
      }

      // Auto-rotate
      autoRotate += 0.002
      torus.rotation.y = autoRotate
      wireframe.rotation.y = autoRotate
      line.rotation.y = autoRotate
      point.rotation.y = autoRotate

      torus.rotation.x = Math.sin(autoRotate * 0.3) * 0.2
      wireframe.rotation.x = Math.sin(autoRotate * 0.3) * 0.2
      line.rotation.x = Math.sin(autoRotate * 0.3) * 0.2
      point.rotation.x = Math.sin(autoRotate * 0.3) * 0.2

      renderer.render(scene, camera)
    }
    animate()

    // Handle window resize
    const handleResize = () => {
      camera.aspect = window.innerWidth / window.innerHeight
      camera.updateProjectionMatrix()
      renderer.setSize(window.innerWidth, window.innerHeight)
    }
    window.addEventListener("resize", handleResize)

    return () => {
      window.removeEventListener("resize", handleResize)
      cancelAnimationFrame(animationId)
      mountRef.current?.removeChild(renderer.domElement)
      renderer.dispose()
    }
  }, [slope, maxWraps, isPlaying, speed])

  const handleReset = () => {
    setIsPlaying(true)
    // Force re-render by changing a dependency
    setMaxWraps(m => m)
  }

  const presetSlopes = [
    { name: "√2", value: Math.sqrt(2) },
    { name: "φ", value: (1 + Math.sqrt(5)) / 2 },
    { name: "π", value: Math.PI },
    { name: "e", value: Math.E },
    { name: "√3", value: Math.sqrt(3) }
  ]

  return (
    <div
      style={{
        position: "relative",
        width: "100%",
        height: "100vh",
        overflow: "hidden"
      }}
    >
      <div ref={mountRef} />
      <div
        style={{
          position: "absolute",
          top: 20,
          left: 20,
          background: "rgba(10, 10, 10, 0.85)",
          color: "white",
          padding: "20px",
          borderRadius: "12px",
          boxShadow: "0 8px 32px rgba(0,0,0,0.5)",
          backdropFilter: "blur(10px)",
          border: "1px solid rgba(255,255,255,0.1)",
          maxWidth: "320px"
        }}
      >
        <h3
          style={{ margin: "0 0 15px 0", fontSize: "22px", fontWeight: "600" }}
        >
          Line Wrapping on Torus
        </h3>

        <div style={{ marginBottom: "20px" }}>
          <div style={{ display: "flex", gap: "10px", marginBottom: "10px" }}>
            <button
              onClick={() => setIsPlaying(!isPlaying)}
              style={{
                flex: 1,
                padding: "10px",
                background: isPlaying ? "#ff5252" : "#4caf50",
                color: "white",
                border: "none",
                borderRadius: "6px",
                cursor: "pointer",
                fontSize: "14px",
                fontWeight: "600"
              }}
            >
              {isPlaying ? "⏸ Pause" : "▶ Play"}
            </button>
            <button
              onClick={handleReset}
              style={{
                flex: 1,
                padding: "10px",
                background: "#2196f3",
                color: "white",
                border: "none",
                borderRadius: "6px",
                cursor: "pointer",
                fontSize: "14px",
                fontWeight: "600"
              }}
            >
              🔄 Reset
            </button>
          </div>
        </div>

        <div style={{ marginBottom: "18px" }}>
          <label
            style={{
              display: "block",
              marginBottom: "8px",
              fontSize: "14px",
              color: "#ccc"
            }}
          >
            Slope:{" "}
            <span style={{ color: "#ff5252", fontWeight: "600" }}>
              {slope.toFixed(6)}
            </span>
          </label>
          <input
            type="range"
            min="1"
            max="5"
            step="0.01"
            value={slope}
            onChange={e => setSlope(parseFloat(e.target.value))}
            style={{ width: "100%", accentColor: "#ff5252" }}
          />
        </div>

        <div style={{ marginBottom: "18px" }}>
          <label
            style={{
              display: "block",
              marginBottom: "8px",
              fontSize: "14px",
              color: "#ccc"
            }}
          >
            Speed:{" "}
            <span style={{ color: "#ffeb3b", fontWeight: "600" }}>
              {speed}x
            </span>
          </label>
          <input
            type="range"
            min="0.5"
            max="5"
            step="0.5"
            value={speed}
            onChange={e => setSpeed(parseFloat(e.target.value))}
            style={{ width: "100%", accentColor: "#ffeb3b" }}
          />
        </div>

        <div style={{ marginBottom: "18px" }}>
          <label
            style={{
              display: "block",
              marginBottom: "8px",
              fontSize: "14px",
              color: "#ccc"
            }}
          >
            Max Wraps:{" "}
            <span style={{ color: "#2196f3", fontWeight: "600" }}>
              {maxWraps}
            </span>
          </label>
          <input
            type="range"
            min="20"
            max="200"
            step="10"
            value={maxWraps}
            onChange={e => setMaxWraps(parseInt(e.target.value))}
            style={{ width: "100%", accentColor: "#2196f3" }}
          />
        </div>

        <div style={{ marginBottom: "10px" }}>
          <div style={{ fontSize: "13px", marginBottom: "8px", color: "#ccc" }}>
            Irrational slopes:
          </div>
          <div style={{ display: "flex", flexWrap: "wrap", gap: "6px" }}>
            {presetSlopes.map(preset => (
              <button
                key={preset.name}
                onClick={() => setSlope(preset.value)}
                style={{
                  padding: "6px 12px",
                  background:
                    Math.abs(slope - preset.value) < 0.001 ? "#ff5252" : "#333",
                  color: "white",
                  border: "1px solid rgba(255,255,255,0.2)",
                  borderRadius: "6px",
                  cursor: "pointer",
                  fontSize: "12px",
                  fontWeight: "600"
                }}
              >
                {preset.name}
              </button>
            ))}
          </div>
        </div>

        <div
          style={{
            fontSize: "11px",
            color: "#999",
            marginTop: "15px",
            lineHeight: "1.5"
          }}
        >
          <p style={{ margin: "5px 0" }}>
            The <span style={{ color: "#ff5252" }}>red line</span> wraps around
            the torus with irrational slope, never closing and densely filling
            the surface.
          </p>
          <p style={{ margin: "5px 0" }}>
            <span style={{ color: "#ffeb3b" }}>Yellow dot</span> shows current
            position.
          </p>
        </div>
      </div>
    </div>
  )
}

export default TorusIrrationalLine
